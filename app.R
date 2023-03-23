library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(purrr)
library(svglite)
library(patchwork)
library(r3dmol)

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
    #tabsetPanel(type = "tabs",
    #            tabPanel("1D Plot", fluid = TRUE, style = "overflow-y:scroll; max-height: 800px",
    #                     br(),
                         box(width = 6, title = "Configure",status = "primary", solidHeader = F, 
                             tabsetPanel(type = "tabs",
                               tabPanel("Inputs", fluid = TRUE,
                                 br(),
                                 textInput("pname","Protein Name",value = "DNMT3A1"),
                                 numericInput("plength","Protein Length",value = 912),
                                 uiOutput("plotting")
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
                             ),
                             br(),
                             actionButton("plot","Plot",
                                          icon = icon(name = "play", lib = "font-awesome",style="color:white"),
                                          style="background-color: #5ABCB9")
                         ),
                         box(width = 6, title = "3D plot",status = "primary", solidHeader = F,
                             #HTML("<script src='http://3Dmol.csb.pitt.edu/build/3Dmol-nojquery.js'></script>"),
                             #HTML("<script src='https://3Dmol.org/build/3Dmol-min.js'></script>"),
                             #HTML("<script src='https://3Dmol.org/build/3Dmol.ui-min.js'></script>"),
                             #HTML('<div style="height: 400px; width: 400px; position: relative;" class="viewer_3Dmoljs" data-href="https://alphafold.ebi.ac.uk/files/AF-Q9Y6K1-F1-model_v4.pdb" data-backgroundcolor="0xffffff" data-select1="resi:9,23,25,28" data-style1="sphere:radius=6,color=blue" data-select2="resi:4,7,10,19" data-style2="sphere:radius=2,color=red"></div>')
                             #htmlOutput("mol")
                             uiOutput("molplotting"),
                             r3dmolOutput("mol")
                         ),
                         box(width = 12, title= "1D plot",status = "primary", solidHeader = F,
                             uiOutput("pplot")
                         )
                #),
                #tabPanel("Undefined", fluid = TRUE, style = "overflow-y:scroll; max-height: 800px",
                         #br()
                         #withSpinner(DT::dataTableOutput("event_summary_table"))
                #)
    #)
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ## Plotting colours
  domain_fill = "#5ABCB9"
  cols<-c("#777DA7","#FE5F55","#99C24D","#FFB703","#219EBC","#1D3461","#885053","#FB8500")
  
  cols_3d = RColorBrewer::brewer.pal(name="Blues",n=9)[4:9]
  
  output$plotting = renderUI({
    tagList(
      numericInput("xmin","X-minimum",0,min = 0),
      numericInput("xmax","X-maximum",input$plength,min = 1),
      numericInput("breaks","X-breaks",50),
      numericInput("vdvp_win","vD/vP window size",value = 0.02),
      helpText("Window size is a fraction of the protein length")
    )
  })
  
  output$molplotting = renderUI({
    tagList(
      textInput("pdbid","PDB ID","Q9Y6K1"),
      checkboxInput("spin","Spin",value = F),
      checkboxInput("labels","Labels",value = F),
      numericInput("first","First Residue",value = 1,max = input$plength,min=1),
      numericInput("last","Last Residue",input$plength,max = input$plength,min=1),
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
      select(Position,Original_Res,New_Res,Allele_Frequency=`Allele Frequency`) |> 
      unique()
  })
  
  ## VdVP calculations
  vdvp = reactive({
  
    index = gnomad() |> mutate(Position=as.numeric(Position)) |>
      complete(Position=1:input$plength) |> pull(Position) |>
      unique()
  
    mut_index = gnomad() |> pull(Position) |> unique()
    
    window = floor(input$plength*input$vdvp_win)
    vp = length(mut_index)/input$plength
    
    vdvp = index %>% map(function(x){
      vd = sum(mut_index %in% x:(x+window)) / window
      vd/vp
    }) |> unlist()
    
    #data <- data.frame(x=index,y=vdvp)
    data <- as.data.frame(spline(index, vdvp))
    data
    
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
                scale_x_continuous(breaks = seq(0,input$plength,input$breaks),labels = seq(0,input$plength,input$breaks))+
                theme(axis.text.x = element_text(vjust = 1,angle = 90,size = 15),
                      legend.text = element_text(size = 20),
                      legend.title = element_text(size = 25),
                      plot.title = element_text(size = 25),
                      axis.title.x=element_blank(),
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
                labs(title = input$pname)
    
    if(nrow(arch())>0){
      p = p + 
        geom_rect(data=arch(),aes(xmin=Start,xmax=End,ymin=-0.5,ymax=0.5),fill=domain_fill) + 
        geom_text(data=arch(),aes(x=Start+((End-Start)/2), y=0, label=Name),size=4,colour="white")
    }
    
    p2 =  vdvp() |> ggplot(aes(x = x, y = y)) +
      theme_classic() +
      geom_line(size = 0.4) +
      scale_x_continuous(breaks=seq(0,input$plength,input$breaks)) +
      labs(x = "Residue", y = "Vd/Vp") +
      theme(axis.text.x = element_text(vjust = 1,angle = 90,size = 15),
            axis.text.y = element_text(hjust = 1,size = 15),
            axis.title.x = element_text(margin = margin(t = 15), size = 20),
            axis.title.y = element_text(margin = margin(r = 15), size = 20))
    
    p/p2+plot_layout(heights = c(3,1))
    
  })  
  
  output$proteinPlot = renderPlot({
    proteinPlot1d() 
  })
  
  ## Download protein plot file
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
  
  
  ## Get selections from Allele Frequency log_scores from Genomad data
  selections = reactive({
    
    s = gnomad() |> select(Position,Allele_Frequency) |>
      filter(Position >= input$first, Position <= input$last) |> ## User defined region
      summarise(Allele_Frequency=sum(Allele_Frequency),.by = Position) |>
      mutate(Log_score=log10(Allele_Frequency*(1e6)),
             Selection = case_when(Log_score <= 1 ~ 1,
                                   Log_score >1 & Log_score <= 2 ~ 2,
                                   Log_score >2 & Log_score <= 3 ~ 3,
                                   Log_score >3 & Log_score <= 4 ~ 4,
                                   Log_score >4 & Log_score <= 5 ~ 5,
                                   Log_score >5 & Log_score <= 6 ~ 6)) ## What about > 6?
    
    list(s1= s |> filter(Selection==1) |> pull(Position),
         s2= s |> filter(Selection==2) |> pull(Position),
         s3= s |> filter(Selection==3) |> pull(Position),
         s3= s |> filter(Selection==4) |> pull(Position),
         s4= s |> filter(Selection==5) |> pull(Position),
         s6= s |> filter(Selection==6) |> pull(Position))
    
  })
  
  ## Get the 3D model from AF database
  r3mol = reactive({
    
    pid = input$pdbid
    
    model=r3dmol(                         # Set up the initial viewer
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 1000
      )
    ) %>%
    m_add_model(                 # Add model to scene
        data = m_bio3d(bio3d::read.pdb(paste0("https://alphafold.ebi.ac.uk/files/AF-",pid,"-F1-model_v4.pdb"))),
        format = "pdb"
    )
    
    model
    
  })
  
  ## Generate 3DMol with R3dmol
  output$mol <- renderR3dmol({
  
    sel = selections()
    pid = input$pdbid
    
    r3d = r3mol() %>%
      m_zoom_to() %>%               # Zoom to encompass the whole scene
      m_add_style(                  
        sel = m_sel(resi = sel[[1]]),      
        style = m_style_sphere(
          color = cols_3d[1],
          colorScheme = "prop",
          radius = 1.15
        )
      ) %>%
      m_set_style(                  
        sel = m_sel(resi = sel[[2]]),      
        style = m_style_sphere(
          color = cols_3d[2],
          colorScheme = "prop",
          radius = 1.5
        )
      ) %>%
      m_set_style(                  
        sel = m_sel(resi = sel[[3]]),      
        style = m_style_sphere(
          color = cols_3d[3],
          colorScheme = "prop",
          radius = 1.85
        )
      ) %>%
      m_set_style(                  
        sel = m_sel(resi = sel[[4]]),      
        style = m_style_sphere(
          color = cols_3d[4],
          colorScheme = "prop",
          radius = 2.15
        )
      ) %>%
      m_set_style(                  
        sel = m_sel(resi = sel[[5]]),      
        style = m_style_sphere(
          color = cols_3d[5],
          colorScheme = "prop",
          radius = 2.5
        )
      ) %>%
      m_set_style(                  
        sel = m_sel(resi = sel[[6]]),      
        style = m_style_sphere(
          color = cols_3d[6],
          colorScheme = "prop",
          radius = 2.85
        )
      ) 
    
      if(input$labels==T){
        r3d = r3d %>%
          m_add_res_labels(
            sel=m_sel(resi = sel[[1]]),
            style = m_style_label(
              backgroundColor = cols_3d[1],
              inFront = T,
              fontSize = 12,
              showBackground = T)
          ) %>%
          m_add_res_labels(
            sel=m_sel(resi = sel[[2]]),
            style = m_style_label(
              backgroundColor = cols_3d[2],
              inFront = T,
              fontSize = 12,
              showBackground = T)
          ) %>%
          m_add_res_labels(
            sel=m_sel(resi = sel[[3]]),
            style = m_style_label(
              backgroundColor = cols_3d[3],
              inFront = T,
              fontSize = 12,
              showBackground = T)
          ) %>%
          m_add_res_labels(
            sel=m_sel(resi = sel[[4]]),
            style = m_style_label(
              backgroundColor = cols_3d[4],
              inFront = T,
              fontSize = 12,
              showBackground = T)
          ) %>%
          m_add_res_labels(
            sel=m_sel(resi = sel[[5]]),
            style = m_style_label(
              backgroundColor = cols_3d[5],
              inFront = T,
              fontSize = 12,
              showBackground = T)
          ) %>%
          m_add_res_labels(
            sel=m_sel(resi = sel[[6]]),
            style = m_style_label(
              backgroundColor = cols_3d[6],
              inFront = T,
              fontSize = 12,
              showBackground = T)
          )
      }

    
      if(input$spin==T){
        r3d = r3d %>% m_spin(speed = 0.3)
      }
      
      r3d

      ## Add labels on and off check box
      ## Min max selection (if required)
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
