library(shiny)
library(bslib)
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
library(forcats)
library(UniprotR)
library(bio3d)
library(ggnewscale)
library(httr2)
library(jsonlite)
library(bslib)


# Define UI
ui<-navbarPage(title="Gnomadic",fluid = T,theme = bs_theme(version = 4, bootswatch = "lumen"), 
               #includeCSS("www/custom.css"),
               tabPanel("1D Plot",fluid = TRUE, style = "overflow-y:scroll; max-height: 1200px",
                        fluidRow(
                          column(width = 3,
                                 h4("Configure"),
                                 br(),
                                 tabsetPanel(type = "tabs",
                                             tabPanel("Inputs", fluid = TRUE,
                                                      br(),
                                                      textInput("pid","Uniprot ID",value = "Q9Y6K1"),
                                                      uiOutput("pinfo"),
                                                      #numericInput("plength","Protein Length",value = 912),
                                                      uiOutput("plotting")
                                             ),
                                             tabPanel("Uploads", fluid = TRUE,
                                                      br(),
                                                      selectInput("gnomad_upload",label = "Gnomad upload",choices = c("Upload","API")),
                                                      uiOutput("gnomad_ui"),
                                                      fileInput(inputId = "consurf_file",label = "Upload Consurf pdb file"),
                                                      uiOutput("fellsUI")
                                             ),
                                             tabPanel("Domains", fluid = TRUE,
                                                      br(),
                                                      selectInput("arch_upload",label = "Protein domain upload",choices = c("Upload","Uniprot")),
                                                      uiOutput("arch_ui"),
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
                                             tabPanel("Downloads", fluid = TRUE,
                                                      br(),
                                                      uiOutput("downloads")
                                             ),
                                 ),
                          ),
                          column(width = 8,
                                 h4("1D Plot"),
                                 tags$div(style="height: 800px",
                                          uiOutput("pplot"),
                                 )
                          ),
                          column(width = 1,
                                 actionButton("plot","Plot",
                                              icon = icon(name = "play", lib = "font-awesome",style="color:white"),
                                              style="background #5ABCB9; background-color: #5ABCB9; border-color: #5ABCB9"),
                                 br(),
                          )
                        )
               ),
               tabPanel("3D Plot",fluid = TRUE, style = "overflow-y:scroll; max-height: 1200px",
                        fluidRow(
                          column(width = 3,
                                 h4("Configure"),
                                 uiOutput("molplotting"),
                                 uiOutput("plot_settings_3d")
                          ),
                          column(width = 7,
                                 h4("3D Plot"),
                                 #HTML("<script src='https://3Dmol.org/build/3Dmol-min.js'></script>"),
                                 #HTML("<script src='https://3Dmol.org/build/3Dmol.ui-min.js'></script>"),
                                 #HTML('<div style="height: 400px; width: 400px; position: relative;" class="viewer_3Dmoljs" data-href="https://alphafold.ebi.ac.uk/files/AF-Q9Y6K1-F1-model_v4.pdb" data-backgroundcolor="0xffffff" data-select1="resi:9,23,25,28" data-style1="sphere:radius=6,color=blue" data-select2="resi:4,7,10,19" data-style2="sphere:radius=2,color=red"></div>')
                                 #htmlOutput("mol")
                                 r3dmolOutput("r3dmol",width = "800px",height = "600px")
                          ),
                          column(width = 2,
                                 downloadButton("pymol", "PyMol Script"),
                                 br(),
                          )
                        )
               )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #bs_themer()
  
  ## Plotting colours
  domain_fill = "#5ABCB9"
    cols<-c("#777DA7","#FE5F55","#99C24D","#FFB703","#219EBC","#1D3461","#885053","#FB8500")
    ms_col="#777DA7"
      consurf_cols = c("1"="#10C8D2","2"="#89FDFD","3"="#D8FDFE","4"="#EAFFFF","5"="#FFFFFF","6"="#FBECF1","7"="#FAC9DE","8"="#F27EAB","9"="#A22664")
      
      ## FELLS Colour
      dh_cols = c(Disorder="#FF0000", HCA="#0F0F0F")
      hsc_cols = c(Helix="#99004D",Strand="#FF9900",Coil="#878787")
      
      ## 3D colours  
      cols_3d = RColorBrewer::brewer.pal(name="Blues",n=9)[4:9]
      
      output$pinfo = renderUI({
        uniprot_info = GetProteinAnnontate(input$pid,columns = c("length","gene_names"))
        tagList(
          textInput("pname","Protein Name",uniprot_info[2]),
          numericInput("plength","Protein Length",uniprot_info[1],min=1)
        )
      })
      
      output$gnomad_ui = renderUI({
        gnomad_trans = NULL
        if(input$gnomad_upload=="API"){
          gnomad_trans = gnomad_api() |> pull(`Transcript ID`) |> unique()
          mane = gnomad_canonical_trans()
          tagList(
            selectInput("transcript",label = "Transcript",choices = gnomad_trans,
                        selected = mane,
                        multiple = F)
          )
        }
        else{
          tagList(
            fileInput(inputId = "gnomad_file",label = "Upload Gnomad csv file")
          )
        }
      })
      
      output$arch_ui = renderUI({
        if(input$arch_upload=="Upload"){
          tagList(
            fileInput(inputId = "arch_file",label = "Upload Protein Domain File")
          )
        }
      })
      
      output$plotting = renderUI({
        tagList(
          numericInput("xmin","X-minimum",0,min = 0),
          numericInput("xmax","X-maximum",input$plength,min = 1),
          numericInput("breaks","X-breaks",50),
          numericInput("vdvp_win","vD/vP Window Size",value = 0.02),
          helpText("Window size: Integer = fixed length, < 1 = fraction of the protein length")
        )
      })
      
      output$molplotting = renderUI({
        tagList(
          textInput("pid2","Uniprot ID","Q9Y6K1"),
          checkboxInput("spin","Spin",value = F),
          checkboxInput("labels","Labels",value = F),
          numericInput("first","First Residue Selection",value = 1,max = input$plength,min=0),
          numericInput("last","Last Residue Selection",input$plength,max = input$plength,min=0),
          actionButton("selectSpheres","Update Selection"),
          br(),
        )    
      })
      
      output$downloads = renderUI({
        tagList(
          textInput("filename_pplot","Plot File Name",value = paste0(input$pid,"_protein_plot.pdf")),
          downloadButton("save_pplot", "Download Plot"),
          helpText("Specify file format as suffix (.pdf,.svg,.png)"),
          textInput("filename_gnomad","Plot File Name",value = paste0(input$pid,"_gnomad.tsv")),
          downloadButton("save_gnomad", "Download Gnomad Data"),
          textInput("filename_vdvp","VdVP File Name",value = paste0(input$pid,"_VdVp.csv")),
          downloadButton("save_vdvp", "Download VdVp"),
          textInput("filename_fells","FELLS File Name",value = paste0(input$pid,"_fells.rds")),
          downloadButton("save_fells", "Download FELLS Data")
        )
      })
      
      ####### FELLS
      
      fells_job = reactiveVal()
      fells_result = reactiveVal()
      fells_class = reactiveVal("btn-primary")
      
      output$fellsUI <- renderUI(
        tagList(
          actionButton("fells",label = "Submit to FELLS",class = fells_class())
        )
      )
      
      observeEvent(input$fells,{
        
        seq=GetProteinAnnontate(input$pid,columns = c("sequence"))
        pseq = paste0(">",input$pid,"\n",seq)
        
        req <- request("http://protein.bio.unipd.it/fellsws/submit") |> 
          req_headers("Content-Type" = "multipart/form-data") |>
          req_body_multipart(sequence = pseq)
        
        resp <- req_perform(req)
        rb = resp |> resp_body_json()
        fells_job(rb$jobid)
        
        showModal(modalDialog(title = "FELLS job submitted"))
        
        job = fells_job()
        status = "submitted"
        fells_class("btn-warning")
        updateActionButton(inputId = "fells",label = "FELLS submitted")
        while(!status=="done"){
          Sys.sleep(10) 
          req <- request(paste0("http://protein.bio.unipd.it/fellsws/status/",job))
          resp = req_perform(req)
          rb = resp |> resp_body_json()
          status = rb$status
        }
        
        if(status=="done"){
          result = rb$names[[1]][[2]]
          
          req2 <- request(paste0("http://protein.bio.unipd.it/fellsws/result/",result))
          res2 = req_perform(req2)
          rb2 = res2 |> resp_body_json()
          fells_result(rb2)
          
          showModal(modalDialog(title = "FELLS job complete"))
          updateActionButton(inputId = "fells",label = "FELLS complete")
          fells_class("btn-success")
        }
        else({
          showModal(modalDialog(title = "FELLS job failed"))
          updateActionButton(inputId = "fells",label = "FELLS failed")
          fells_class("btn-danger")
        })
        
      })
      
      fells_hd <- reactive({
        
        fr = fells_result()
        
        hca = fr$hca |> unlist()
        dis = fr$p_dis |> unlist()
        
        df <- data.frame(index = seq(1,input$plength,1),
                         HCA=as.numeric(hca),
                         Disorder=as.numeric(dis)) |>
          mutate(Disorder = Disorder * -1) |>
          pivot_longer(-index) |> 
          mutate(alpha=case_when(name=="Disorder"~value*-1,.default=value)) |>
          mutate(alphag=cut_interval(alpha,n=5))
        
        df
      })
      
      fells_hsc <- reactive({
        
        fr = fells_result()
        
        phelix = fr$p_h |> unlist()
        pstrand = fr$p_e |> unlist()
        pcoil = fr$p_c |> unlist()
        
        df <- data.frame(index = seq(1,input$plength,1),
                         Helix=as.numeric(phelix),
                         Strand=as.numeric(pstrand),
                         Coil=as.numeric(pcoil)) |>
          mutate(Coil = Coil * -1) |>
          pivot_longer(-index) |> 
          mutate(alpha=case_when(name=="Coil"~value*-1,.default=value)) |>
          mutate(alphag=cut_interval(alpha,n=5))
        df
      })
      
      #############
      
      ## Extra settings for 3dmol
      output$plot_settings_3d <- renderUI({
        tagList(
          br(),
          checkboxInput("surface","Add surface"),
          selectInput(
            inputId = "set_style",
            label = "Set main style",
            choices = c("Line", "Cross", "Stick", "Sphere", "Cartoon"),
            selected = "Line"
          ),
          sliderInput(
            inputId = "set_slab",
            label = "Set slab of view",
            min = -150,
            value = c(-150, 150),
            animate = TRUE,
            step = 10,
            max = 150,
            dragRange = TRUE
          )#,
          #radioButtons(
          #  label = "Set view projection scheme",
          #  inputId = "set_projection",
          #  choices = c("perspective", "orthographic"),
          #  inline = TRUE
          #),
          #sliderInput(
          #  inputId = "set_perceived_distance",
          #  label = "Set perceived distance",
          #  min = 0,
          #  max = 500,
          #  value = 300
          #)
        )
      })
      
      observeEvent(input$plot,{
        output$pplot = renderUI({
          tagList(
            withSpinner(plotOutput("proteinPlot"))
          )
        })
      })
      
      ## Amino acid code
      pct = read_tsv("protein_code.tsv",col_names = T)
      pc = pct$One
      names(pc) = pct$Three
      
      ## Hard coded annotation        
      #arch2 = read_tsv("test_data/Inputs/4_DNMT3A1_architecture_file.txt")
      #anno2 = read_tsv("test_data/Inputs/4_DNMT3A1_post_translation_file.txt") |> 
      #  mutate(Annotation="PTM")
      
      ## Gnomad Filters        
      clinVar = c("Uncertain significance","Likely pathogenic","Pathogenic",
                  "Conflicting interpretations of pathogenicity")
      
      ## Gnomad API query
      query <- reactive({
        paste0('{\n  gene(gene_symbol: "', input$pname, '", reference_genome: GRCh38) {
                  symbol
                  name
                  canonical_transcript_id
                  variants(dataset: ',"gnomad_r4",') {
                    variant_id
                    transcript_id
                    chrom
                    pos
                    ref
                    alt
                    consequence
                    hgvsp
                    genome{
                      ac
                      an
                    }
                    exome{
                      ac
                      an
                    }
                  }
                  clinvar_variants {
                    variant_id
                    clinical_significance
                  }
                }
          }'
        )
      })
      
      gnomad_canonical_trans <- reactiveVal()
      
      gnomad_api <- reactive({
        jsondata <- request("https://gnomad.broadinstitute.org/api/?") %>%
          req_body_json(list(query=query(), variables="null")) %>%
          req_perform() %>%
          resp_body_json()
        
        gnomad_v<-jsondata$data$gene$variants |> map(~unlist(.x) |> t() |> as.data.frame()) |> bind_rows()
        gnomad_c<-jsondata$data$gene$clinvar_variants |> map(~unlist(.x) |> t() |> as.data.frame()) |> bind_rows()
        
        gnomad_canonical_trans(jsondata$data$gene$canonical_transcript_id)
        
        gnomadF <- gnomad_v |> left_join(gnomad_c,by="variant_id") |>
          mutate(across(starts_with(c("genome","exome")),~coalesce(as.numeric(.x),0))) |>
          mutate(af = (genome.ac+exome.ac) / (genome.an+exome.an)) |>
          select(`VEP Annotation` = consequence,
                 `Protein Consequence` = hgvsp,
                 `ClinVar Clinical Significance` = clinical_significance,
                 `Allele Frequency` = af,
                 `Transcript ID` = transcript_id)
        
      }) 
      
      ## Read in gnomad file and filter        
      gnomad = reactive({
        
        gnomadF = NULL
        if(!is.null(input$gnomad_file)){
          tryCatch(
            {
              gnomadF = read_csv(input$gnomad_file$datapath,col_names = T)
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
        else{
          if(input$gnomad_upload=="API"){
            gnomadF = gnomad_api() |> filter(`Transcript ID` == input$transcript)
          }
          else{
            ## Default file
            gnomadF = read_csv("test_data/Inputs/1-2-3-5_DNMT3A1_gnomad_file.csv",col_names = T)
          }
        }
        
        ## Auto-detect v2 and v3 files?
        
        gnomadF |> filter(`VEP Annotation` == "missense_variant",
                          #`Filters - exomes` == "PASS",
                          !`ClinVar Clinical Significance` %in% clinVar) |>
          mutate(`Protein Consequence` = str_remove(`Protein Consequence`,"^p.")) |>
          mutate(Original_Res = pc[str_sub(`Protein Consequence`,1,3)],
                 New_Res = pc[str_sub(`Protein Consequence`,-3)],
                 Position = str_remove_all(`Protein Consequence`,"[A-Z,a-z]") |> as.numeric()) |>
          select(Position,Original_Res,New_Res,Allele_Frequency=`Allele Frequency`) |> 
          unique()
      })
      
      ## VdVP calculations
      vdvp = reactive({
        
        index = 1:input$plength
        mut_index = gnomad() |> pull(Position) |> unique()
        
        if(input$vdvp_win<1){
          window = floor(input$plength*input$vdvp_win) ## Fraction of length
        }
        else{
          window = input$vdvp_win ## Fixed size
        }
        vp = length(mut_index)/input$plength
        
        vdvp = index %>% map(function(x){
          vd = sum(mut_index %in% x:(x+window)) / window
          vd/vp
        }) |> unlist()
        
        ## vd/vp reverse
        #index.rev = input$plength:1
        
        #vdvp2 = index.rev %>% map(function(x){
        #  vd = sum(mut_index %in% x:(x-window)) / window
        #  vd/vp
        #}) |> unlist()
        
        #data <- data.frame(x=index,y=vdvp)
        data <- as.data.frame(spline(index+(window/2),vdvp)) #|> mutate(direction="forward") ## Centre window
        data
        #data2 <- as.data.frame(spline(index.rev-(window/2),vdvp2)) |> mutate(direction="reverse") ## Centre window
        #bind_rows(data,data2)
        
      })
      
      ## Get protein architecture
      arch = reactive({
        
        archF = NULL
        if(!is.null(input$arch_file) & input$arch_upload == "Upload"){
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
        else{
          if(input$arch_upload == "Uniprot"){
            fd <- GetFamily_Domains(input$pid)
            if(!is.na(fd$Domain..FT.)){
              archF <- fd$Domain..FT. |> str_split_1(pattern = "DOMAIN") |> map(~str_split_1(.x,";")) |> 
                map(function(x){
                  if(x[1]!=""){
                    coord=x[1] |> str_remove_all(" ") |> str_split_fixed("\\.\\.",2)
                    Start=coord[1]
                    End=coord[2]
                    Name=x[2] |> str_remove_all(" ") |> str_remove("/note=")
                    data.frame(Name,Start,End) |>
                      mutate(Start = as.numeric(Start),
                             End = as.numeric(End))
                  }
                }) |>
                bind_rows()
            }
          }
        }
        
        archT = paste0(input$arch,"\n")
        if(!archT == "\n"){
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
          mutate(Annotation=as.factor(Annotation) %>% fct_relevel("MissenseVariant")) |>
          group_by(Annotation) |>
          mutate(y = cur_group_id()-1)
        
      })
      
      consurf = reactive({
        
        pdbF = NULL
        if(!is.null(input$consurf_file)){
          tryCatch(
            {
              pdbF = read.pdb(input$consurf_file$datapath,ATOM.only = T)
              pdbF$atom |> select(resid,resno,score=b) |> unique()
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
        else{
          NULL
          ## Default file
          #pdbF = read.pdb("test_data/AF-Q12906-F1-model_v4_ATOMS_section_With_ConSurf.pdb",ATOM.only = T)
        }
        
      })
      
      
      ## Plot protein
      proteinPlot1d = eventReactive(input$plot,{
        
        p = annotation() |> ggplot() +
          geom_rect(xmin = 0,xmax=input$plength,ymin=-0.1,ymax=0.1,fill="darkgrey")+
          theme_bw()+
          scale_x_continuous(breaks = seq(0,input$plength,input$breaks),labels = seq(0,input$plength,input$breaks))+
          theme(axis.text.x = element_text(vjust = 1,angle = 90,size = 15),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 20),
                plot.title = element_text(size = 25),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.margin = margin(2, 2, 2, 2, "cm")
          ) +
          geom_point(aes(x=Start,y=y+1,colour=Annotation),alpha=0.8,shape="|",size=7)+
          #geom_tile(aes(x=Start,y=y+1.5,colour=Annotation,fill=Annotation),alpha=0.5)+
          scale_colour_manual(values = cols)+
          coord_cartesian(ylim=c(-1,max(annotation()$y)+2),xlim=c(input$xmin,input$xmax))+
          labs(title = input$pname)
        
        if(!is.null(consurf())){
          p = p + 
            new_scale_colour()+
            geom_point(data = consurf(), aes(x=resno,colour=score,y=-1),alpha=0.8,shape="|",size=7)+
            scale_colour_gradientn(colours = consurf_cols,guide = guide_colourbar(direction = "horizontal",title.position = "top",label.theme = element_text(size = 10))) +
            labs(colour = "Consurf")+
            #scale_colour_manual(values = consurf_cols,)+
            coord_cartesian(ylim=c(-2,max(annotation()$y)+2),xlim=c(input$xmin,input$xmax))
        }
        
        if(nrow(arch())>0){
          p = p + 
            geom_rect(data=arch(),aes(xmin=Start,xmax=End,ymin=-0.5,ymax=0.5),fill=domain_fill) + 
            geom_text(data=arch(),aes(x=Start+((End-Start)/2), y=0, label=Name),size=4,colour="white")
        }
        
        p2 =  vdvp() |> ggplot(aes(x = x, y = y)) +  #group=direction,colour=direction)) +
          theme_classic() +
          geom_line(linewidth = 0.4,colour="#19647E") +
          #scale_colour_manual(values = c("#19647E","#749C75")) +
          scale_x_continuous(breaks=seq(0,input$plength,input$breaks)) +
          labs(x = "Residue", y = "Vd/Vp") +
          theme(axis.text.x = element_text(vjust = 1,angle = 90,size = 15),
                axis.text.y = element_text(hjust = 1,size = 15),
                axis.title.x = element_text(margin = margin(t = 15), size = 20),
                axis.title.y = element_text(margin = margin(r = 15), size = 20))+
          coord_cartesian(xlim=c(input$xmin,input$xmax))
        
        if(fells_class()=="btn-success"){
          p3 <- fells_hsc() |> ggplot(aes(index,value,fill=name,alpha=alpha))+
            geom_col(position = "identity") +
            scale_fill_manual(values = hsc_cols)+
            guides(alpha="none")+
            theme_minimal()
          
          p4 <- fells_hd() |> ggplot(aes(index,value,fill=name,alpha=alpha))+
            geom_col(position = "identity") +
            scale_fill_manual(values = dh_cols)+
            guides(alpha="none")+
            theme_minimal()
          
          p/p2/p3/p4+plot_layout(heights = c(4,1,1,1))
        }
        else{
          p/p2+plot_layout(heights = c(3,1))
        }
      })  
      
      output$proteinPlot = renderPlot({
        proteinPlot1d() 
      },height=700)
      
      ## Download protein plot file
      output$save_pplot <- downloadHandler(
        filename = function() {
          input$filename_pplot
        },
        content = function(file) {
          if(!is.null(proteinPlot1d())){
            p <- proteinPlot1d()
            save_height = 5+max(annotation()["y"])
            if(fells_class()=="btn-success"){
              save_height = save_height + 10
            }
            ggsave(filename = file,plot = p,width = 20, height = save_height)
          }
        }
      )
      
      ## Download Gnomad file
      output$save_gnomad <- downloadHandler(
        filename = function() {
          input$filename_gnomad
        },
        content = function(file) {
          if(!is.null(proteinPlot1d())){
            gf <- gnomad()
            write_tsv(gf, file = file)
          }
        }
      )
      
      ## Download FELLS file
      output$save_fells <- downloadHandler(
        filename = function() {
          input$filename_fells
        },
        content = function(file) {
          if(fells_class()=="btn-success"){
            fr <- fells_result()
            saveRDS(fr,file = file)
          }
        }
      )
      
      ## Download vdvp file
      output$save_vdvp <- downloadHandler(
        filename = function() {
          input$filename_vdvp
        },
        content = function(file) {
          if(!is.null(proteinPlot1d())){
            v <- vdvp()
            write_csv(v, file = file)
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
      
      ## Get the 3D model from AF database and generate model
      output$r3dmol = renderR3dmol({
        
        if(!is.null(input$pid2)){
          pid = input$pid2
          
          r3dmol(                         # Set up the initial viewer
            viewer_spec = m_viewer_spec(
              cartoonQuality = 25,
              lowerZoomLimit = 5,
              upperZoomLimit = 1000
            )
          ) %>%
            m_add_model(                 # Add model to scene
              data = m_bio3d(bio3d::read.pdb(paste0("https://alphafold.ebi.ac.uk/files/AF-",pid,"-F1-model_v4.pdb"))),
              format = "pdb"
            ) %>% 
            m_zoom_to()
        }
      })
      
      observeEvent(input$selectSpheres,{
        
        if(!is.null(input$pid2)){
          
          style <- switch(
            input$set_style,
            "Line" = list(line = list()),
            "Cartoon" = list(cartoon = list()),
            "Stick" = list(stick = list()),
            "Cross" = list(cross = list()),
            "Sphere" = list(sphere = list())
          )
          m_remove_all_shapes(id="r3dmol")
          m_set_style(
            id="r3dmol",
            style = style
          )
          
          if(input$last>0){
            sel = selections()
            
            m_add_style(
              id="r3dmol",
              sel = m_sel(resi = sel[[1]],atom = "CA"),      
              style = m_style_sphere(
                color = cols_3d[1],
                colorScheme = "prop",
                radius = 1.15
              )
            )
            m_add_style(
              id="r3dmol",
              sel = m_sel(resi = sel[[2]],atom = "CA"),      
              style = m_style_sphere(
                color = cols_3d[2],
                colorScheme = "prop",
                radius = 1.5
              )
            )
            m_add_style(
              id="r3dmol",
              sel = m_sel(resi = sel[[3]],atom = "CA"),      
              style = m_style_sphere(
                color = cols_3d[3],
                colorScheme = "prop",
                radius = 1.85
              )
            ) 
            m_add_style(
              id="r3dmol",
              sel = m_sel(resi = sel[[4]],atom = "CA"),      
              style = m_style_sphere(
                color = cols_3d[4],
                colorScheme = "prop",
                radius = 2.15
              )
            ) 
            m_add_style(
              id="r3dmol",
              sel = m_sel(resi = sel[[5]],atom = "CA"),      
              style = m_style_sphere(
                color = cols_3d[5],
                colorScheme = "prop",
                radius = 2.5
              )
            ) 
            m_add_style(
              id="r3dmol",
              sel = m_sel(resi = sel[[6]],atom = "CA"),      
              style = m_style_sphere(
                color = cols_3d[6],
                colorScheme = "prop",
                radius = 2.85
              )
            )
          }
        }
      })
      
      ## Turn on/off labels for selections
      observeEvent(input$labels,{
        
        if(!is.null(input$pid2)){
          if(input$labels==T){
            sel = selections()
            
            m_add_res_labels(
              id="r3dmol",
              sel=m_sel(resi = sel[[1]]),
              style = m_style_label(
                backgroundColor = cols_3d[1],
                inFront = T,
                fontSize = 12,
                showBackground = T)
            )
            m_add_res_labels(
              id="r3dmol",
              sel=m_sel(resi = sel[[2]]),
              style = m_style_label(
                backgroundColor = cols_3d[2],
                inFront = T,
                fontSize = 12,
                showBackground = T)
            )
            m_add_res_labels(
              id="r3dmol",
              sel=m_sel(resi = sel[[3]]),
              style = m_style_label(
                backgroundColor = cols_3d[3],
                inFront = T,
                fontSize = 12,
                showBackground = T)
            )
            m_add_res_labels(
              id="r3dmol",
              sel=m_sel(resi = sel[[4]]),
              style = m_style_label(
                backgroundColor = cols_3d[4],
                inFront = T,
                fontSize = 12,
                showBackground = T)
            )
            m_add_res_labels(
              id="r3dmol",
              sel=m_sel(resi = sel[[5]]),
              style = m_style_label(
                backgroundColor = cols_3d[5],
                inFront = T,
                fontSize = 12,
                showBackground = T)
            )
            m_add_res_labels(
              id="r3dmol",
              sel=m_sel(resi = sel[[6]]),
              style = m_style_label(
                backgroundColor = cols_3d[6],
                inFront = T,
                fontSize = 12,
                showBackground = T)
            )
          }
          else{
            m_remove_all_labels(id="r3dmol")
          }
        }
      })
      
      ## Turn on / off spin
      observeEvent(input$spin,{
        if(input$spin==T){
          m_spin(id = "r3dmol",speed = 0.3)
        }
        else{
          m_spin(id = "r3dmol",speed = 0)
        }
      })
      
      observeEvent(input$set_style, {
        style <- switch(
          input$set_style,
          "Line" = list(line = list()),
          "Cartoon" = list(cartoon = list()),
          "Stick" = list(stick = list()),
          "Cross" = list(cross = list()),
          "Sphere" = list(sphere = list())
        )
        
        #sel = selections()
        
        m_set_style(id = "r3dmol",
                    #sel = m_sel(resi = unlist(sel,use.names = F),atom="CA",invert = T),
                    #sel = m_sel(resi = c(1:input$plength)[-c(unlist(sel))],invert = F),
                    style = style
        )
      })
      
      #observeEvent(input$set_projection, {
      #  m_set_projection(id = "r3dmol", scheme = input$set_projection)
      #})
      
      observeEvent(input$set_slab, {
        m_set_slab(id = "r3dmol",
                   near = input$set_slab[1],
                   far = input$set_slab[2])
      })
      
      observeEvent(input$surface,{
        if(input$surface==T){
          m_add_surface(id = "r3dmol",style = m_style_surface(opacity = 0.4))
        }
        else{
          m_remove_all_surfaces(id = "r3dmol")
        }
      })
      
      
      #observeEvent(input$set_perceived_distance, {
      #  m_set_preceived_distance(id = "r3dmol", dist = input$set_perceived_distance)
      #})
      
      ## Generate pyMol script
      pymolScript = reactive({
        sel = selections()
        i = 1:6  
        ps = "bg_color white\ncolor white, DNMT3A1\nutil.performance(0)\nspace rgb\nset ray_shadows,off\n"
        pss = i |> map(~paste0("select s",.x," (((i;",paste(sort(sel[[.x]]),collapse = ","),") and n; CA) and ",input$pid2,"\n"))  |> reduce(paste0)
        ps = paste0(ps,pss,"show spheres, (s1,s2,s3,s4,s5,s6)\nset sphere_scale, 1.15, s1\nset sphere_scale, 1.50, s2\nset sphere_scale, 1.85, s3\nset sphere_scale, 2.15, s4\nset sphere_scale, 2.50, s5\nset sphere_scale, 2.85, s6\nset_color b1, [100,120,250]\nset_color b2, [35,40,200]\nset_color b3, [18,0,150]\nset_color b4, [9,0,100]\nset_color b5, [0,0,50]\nset_color b6, [0,0,0]\ncolor b1, s1\ncolor b2, s2\ncolor b3, s3\ncolor b4, s4\ncolor b5, s5\ncolor b6, s6\n")
        ps
      })
      
      ## Download pymol file
      output$pymol <- downloadHandler(
        filename = function() {
          paste0(input$pid2,".pymol")
        },
        content = function(file) {
          if(!is.null(input$pid2)){
            ps=pymolScript()
            write(ps,file = file)
          }
        }
      )
      
}

# Run the application 
shinyApp(ui = ui, server = server)
