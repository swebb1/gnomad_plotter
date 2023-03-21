library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)

# Define UI
ui<-dashboardPage(
  
  # Application title
  dashboardHeader(title = "GnomadPlot", titleWidth = 240),
  dashboardSidebar(
    br(),
    div(style="text-align: center;",
        #actionButton("load","Load from Spond",icon = icon(name = "upload", lib = "font-awesome",style="color:white"),style="background-color: #367ea9;color:white"),
        br(),
        br()
        #img(src="GlasgowUltimateLogo2013.png",width="80%")
    )
  ),
  dashboardBody(
    includeCSS("www/custom.css"),
    tabsetPanel(type = "tabs",
                tabPanel("1D Plot", fluid = TRUE, style = "overflow-y:scroll; max-height: 800px",
                         br(),
                         box(width = 12, title = "Config",status = "primary", solidHeader = F, 
                             textInput("pname","Protein Name",value = "DNMT3A1"),
                             numericInput("plength","Protein Length",value = 912),
                             textAreaInput("arch","Protein Domains"),
                             helpText("Protein domain: name,start,end separated by comma"),
                             tableOutput("pp"),
                             actionButton('addAnno', 'Add annotations'),
                             div(id = 'placeholder'),
                             br(),
                             actionButton("plot","Plot",
                                          icon = icon(name = "play", lib = "font-awesome",style="color:white"),
                                          style="background-color: #5ABCB9")
                         ),
                         box(width = 12, title= "1D plot",status = "primary", solidHeader = F,
                             withSpinner(plotOutput("proteinPlot")))
                ),
                tabPanel("Summary", fluid = TRUE, style = "overflow-y:scroll; max-height: 800px",
                         br()
                         #withSpinner(DT::dataTableOutput("event_summary_table"))
                )
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ## Plot colours
  domain_fill = "#5ABCB9"
  cols<-c("MissenseVariant"="#777DA7","PTM"="#FE5F55")
            
  ## Amino acid code
  pct = read_tsv("protein_code.tsv",col_names = T)
  pc = pct$One
  names(pc) = pct$Three
  
  ## Hard coded annotation        
  arch2 = read_tsv("../test_data/Inputs/4_DNMT3A1_architecture_file.txt")
  ptm2 = read_tsv("../test_data/Inputs/4_DNMT3A1_post_translation_file.txt") |> 
    mutate(Annotation="PTM")
  
  ## Gnomad Filters        
  clinVar = c("Uncertain significance","Likely pathogenic","Pathogenic",
                      "Conflicting interpretations of pathogenicity")
  
  ## Get protein architecture
  arch = eventReactive(input$plot,{
    archT = input$arch
    read_csv(archT,col_names = c("Name","Start","End"))
  })
  
  ## Print table for user verification
  output$pp = renderTable(arch())
  
  ## Dynamically add annotations
  annoCount <- reactiveVal(0)

  observeEvent(input$addAnno,{

    #btn=sum(input$addAnno, 1)
    annoCount(annoCount()+1)
    btn=annoCount()

    insertUI(
      selector = '#placeholder',
      ui = div(
        id = paste0("div_",btn),
        h5(paste0("Annotation ",btn)),
        textInput(paste0("annoName_",btn),label = "Name"),
        textAreaInput(paste0("annoPos_",btn),label = "Positions"),
        actionButton(paste0("remove_",btn),label = "Remove")
      )
    )
    
    observeEvent(input[[paste0("remove_",btn)]], {
        rmv = paste0("#div_",btn)
        removeUI(
          selector = rmv
        )
        #annoCount(annoCount()-1)
    })
    
  })
  

  ## Read in gnomad file and filter        
  gnomad = reactive({
          read_csv("../test_data/Inputs/1-2-3-5_DNMT3A1_gnomad_file.csv",col_names = T) |>
          select(Filters_exomes = `Filters - exomes`,
                   Protein_Consequence = `Protein Consequence`,
                   VEP_Annotation = `VEP Annotation`,
                   ClinVar_Clinical_Significance = `ClinVar Clinical Significance`) |>
          filter(VEP_Annotation == "missense_variant",
                   Filters_exomes == "PASS",
                   !ClinVar_Clinical_Significance %in% clinVar) |>
          select(-c(VEP_Annotation,Filters_exomes,ClinVar_Clinical_Significance)) |>
          mutate(Protein_Consequence = str_remove(Protein_Consequence,"^p.")) |>
          mutate(Original_Res = pc[str_sub(Protein_Consequence,1,3)],
                   New_Res = pc[str_sub(Protein_Consequence,-3)],
                   Position = str_remove_all(Protein_Consequence,"[A-Z,a-z]")) |>
          select(Position,Original_Res,New_Res)
  })
  
  ## Create table of annotations
  annotation = reactive({
    ## Get missense mutations as annotation
    mut = tibble(site=as.numeric(gnomad()$Position),Annotation="MissenseVariant")
    
    ## Bind missense and other annotations
    bind_rows(mut,ptm2) |>
        group_by(Annotation) %>%
        mutate(y = cur_group_id())
  })
  
  ## Plot protein
  observeEvent(input$plot,{
    output$proteinPlot<-renderPlot({
      arch2 |> ggplot() +
                geom_rect(xmin = 0,xmax=input$plength,ymin=-0.2,ymax=0.2,fill="darkgrey")+
                geom_rect(aes(xmin=Start,xmax=End,ymin=-0.5,ymax=0.5),fill=domain_fill)+
                geom_text(aes(x=Start+((End-Start)/2), y=0, label=Name),size=4,colour="white") +
                theme_bw()+
                scale_x_continuous(breaks = seq(0,input$plength,25),labels = seq(0,input$plength,25))+
                theme(axis.text.x = element_text(vjust = 1,angle = 90),
                      axis.title.x = element_text(margin = margin(t = 10)),
                      axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank()
                ) +
                geom_point(data = annotation(),mapping=aes(x=site,colour=Annotation,y=y+1),alpha=0.5,size=3)+
                scale_colour_manual(values = cols)+
                coord_cartesian(ylim=c(-2,max(annotation()$y)+2),xlim=c(0,input$plength))+
                labs(x="Position",title = input$pname)
    })
  })  
}

# Run the application 
shinyApp(ui = ui, server = server)
