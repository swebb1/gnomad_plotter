library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(purrr)
library(svglite)

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
                         box(width = 6, title = "Configure",status = "primary", solidHeader = F, 
                             tabsetPanel(type = "tabs",
                               tabPanel("Inputs", fluid = TRUE,
                                 br(),
                                 textInput("pname","Protein Name",value = "DNMT3A1"),
                                 numericInput("plength","Protein Length",value = 912)
                               ),
                               tabPanel("Protein Domains", fluid = TRUE,
                                 br(),
                                 fileInput(inputId = "arch_file",label = "Upload Protein Domain File"),
                                 textAreaInput("arch","Paste Protein Domains"),
                                 helpText("Protein domain: name,start,end separated by comma"),
                                 tableOutput("arch_table")
                               ),
                               tabPanel("Annotations", fluid = TRUE,
                                 br(),
                                 fileInput(inputId = "anno_file",label = "Upload Annotation File"),
                                 textAreaInput("anno","Paste Annotations"),
                                 helpText("Annotations: name,position separated by comma"),
                               ),
                               tabPanel("Plotting Arguments", fluid = TRUE,
                                 br(),
                                 uiOutput("plotting")
                               )
                             ),
                             br(),
                             actionButton("plot","Plot",
                                          icon = icon(name = "play", lib = "font-awesome",style="color:white"),
                                          style="background-color: #5ABCB9")
                         ),
                         box(width = 6, title = "3D plot",status = "primary", solidHeader = F,
                             HTML('<script src="https://3Dmol.org/build/3Dmol-min.js"></script>'),
                             HTML('<script src="https://3Dmol.org/build/3Dmol.ui-min.js"></script>'),
                             HTML('<div style="height: 400px; width: 400px; position: relative;" class="viewer_3Dmoljs" data-pdb="6W8B" data-backgroundcolor="0xffffff" data-style="sphere" data-ui="true"></div>')
                         ),
                         box(width = 12, title= "1D plot",status = "primary", solidHeader = F,
                                 uiOutput("pplot")
                         )
                ),
                tabPanel("Undefined", fluid = TRUE, style = "overflow-y:scroll; max-height: 800px",
                         br()
                         #withSpinner(DT::dataTableOutput("event_summary_table"))
                )
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ## Plotting colours
  domain_fill = "#5ABCB9"
  cols<-c("#777DA7","#FE5F55","#99C24D","#FFB703","#219EBC","#1D3461","#885053","#FB8500")
  
  output$plotting = renderUI({
    tagList(
      numericInput("xmin","X-minimum",0,min = 0),
      numericInput("xmax","X-maximum",input$plength,min = 0)
    ) 
  })
  
  observeEvent(input$plot,{
    output$pplot = renderUI({
      tagList(
        withSpinner(plotOutput("proteinPlot")),
        h4("Download"),
        textInput("filename_pplot","File name",value = "Protein_plot.pdf"),
        downloadButton("save_pplot", "Download"),
        helpText("Specify file format as suffix (.pdf,.svg,.png)")
      )
    })
  })
  
  ## Amino acid code
  pct = read_tsv("protein_code.tsv",col_names = T)
  pc = pct$One
  names(pc) = pct$Three
  
  ## Hard coded annotation        
  arch2 = read_tsv("test_data/Inputs/4_DNMT3A1_architecture_file.txt")
  anno2 = read_tsv("test_data/Inputs/4_DNMT3A1_post_translation_file.txt") |> 
    mutate(Annotation="PTM")
  
  ## Gnomad Filters        
  clinVar = c("Uncertain significance","Likely pathogenic","Pathogenic",
                      "Conflicting interpretations of pathogenicity")
  
  
  ## Read in gnomad file and filter        
  gnomad = reactive({
    read_csv("test_data/Inputs/1-2-3-5_DNMT3A1_gnomad_file.csv",col_names = T) |>
      filter(`VEP Annotation` == "missense_variant",
             `Filters - exomes` == "PASS",
             !`ClinVar Clinical Significance` %in% clinVar) |>
      mutate(`Protein Consequence` = str_remove(`Protein Consequence`,"^p.")) |>
      mutate(Original_Res = pc[str_sub(`Protein Consequence`,1,3)],
             New_Res = pc[str_sub(`Protein Consequence`,-3)],
             Position = str_remove_all(`Protein Consequence`,"[A-Z,a-z]")) |>
      select(Position,Original_Res,New_Res)
  })
  
  
  ## Get protein architecture
  arch = reactive({
    
    archF = NULL
    if(!is.null(input$arch_file)){
      tryCatch(
        {
          archF = read_csv(input$arch_file$datapath,col_names = c("Name","Start","End")) |>
            arrange(Start)
        },
        error = function(e){
          message("Cannot upload file")
          showModal(modalDialog(
            title = "Error",
            "Cannot upload file"
          ))
          return(NULL)
        }
      )
    }
    
    archT = input$arch
    if(!archT == ""){
      archT = read_csv(archT,col_names = c("Name","Start","End")) |>
        arrange(Start)
    }
    else{
      archT = NULL
    }
    
    bind_rows(archF,archT)
    
  })
  
  ## Print table for user verification
  output$arch_table = renderTable(arch())
  
  ## Get extra annotations
  anno = reactive({
      
      annoF = NULL
      if(!is.null(input$anno_file)){
        tryCatch(
          {
            annoF = read_csv(input$anno_file$datapath,col_names = c("Annotation","Start")) |>
              arrange(Annotation)
          },
          error = function(e){
            message("Cannot upload file")
            showModal(modalDialog(
              title = "Error",
              "Cannot upload file"
            ))
            return(NULL)
          }
        )
      }
      
      annoT = input$anno
      if(!annoT==""){
        annoT = read_csv(annoT,col_names = c("Annotation","Start")) |>
          arrange(Annotation)
      }
      else{
        annoT = NULL
      }
      
      bind_rows(annoF,annoT)
      
  })
  
  ## Create table of annotations
  annotation = reactive({
    ## Get missense mutations as annotation
    mut = tibble(Annotation="MissenseVariant",Start=as.numeric(gnomad()$Position))
    
    bind_rows(mut,anno()) |>
      group_by(Annotation) |>
      mutate(y = cur_group_id())
    
  })
  

  ## Plot protein
  proteinPlot1d = eventReactive(input$plot,{
      
    p = annotation() |> ggplot() +
                geom_rect(xmin = 0,xmax=input$plength,ymin=-0.2,ymax=0.2,fill="darkgrey")+
                theme_bw()+
                scale_x_continuous(breaks = seq(0,input$plength,25),labels = seq(0,input$plength,25))+
                theme(axis.text.x = element_text(vjust = 1,angle = 90,size = 15),
                      axis.title.x = element_text(margin = margin(t = 15), size = 20),
                      legend.text = element_text(size = 20),
                      legend.title = element_text(size = 25),
                      plot.title = element_text(size = 25),
                      axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank()
                ) +
                geom_point(aes(x=Start,colour=Annotation,y=y+1),alpha=0.5,size=3)+
                scale_colour_manual(values = cols)+
                coord_cartesian(ylim=c(-1,max(annotation()$y)+2),xlim=c(input$xmin,input$xmax))+
                labs(x="Position",title = input$pname)
    
    if(nrow(arch())>0){
      p = p + 
        geom_rect(data=arch(),aes(xmin=Start,xmax=End,ymin=-0.5,ymax=0.5),fill=domain_fill) + 
        geom_text(data=arch(),aes(x=Start+((End-Start)/2), y=0, label=Name),size=4,colour="white")
    }
    
    p
    
  })  
  
  output$proteinPlot = renderPlot({
    proteinPlot1d() 
  })
  
  # Download protein plot file
  output$save_pplot <- downloadHandler(
    filename = function() {
      input$filename_pplot
    },
    content = function(file) {
      if(!is.null(proteinPlot1d())){
        p <- proteinPlot1d()
        ggsave(filename = file,plot = p,width = 20, height = 5+max(annotation()["y"]))
      }
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)
