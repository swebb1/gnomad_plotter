library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)


## Plot colours
domain_fill = "#5ABCB9"
cols<-c("MissenseVariant"="#777DA7","PTM"="#FE5F55")
ptm_fill = "#FE5F55"
mut_col = "#777DA7"
col3 = "#1D3461"
col4 = "#885053"

pct = read_tsv("protein_code.tsv",col_names = T)
pc = pct$One
names(pc) = pct$Three

arch = read_tsv("../test_data/Inputs/4_DNMT3A1_architecture_file.txt")
ptm = read_tsv("../test_data/Inputs/4_DNMT3A1_post_translation_file.txt")
plength = 912
pname = "DNMT3A1"

clinVar = c("Uncertain significance","Likely pathogenic","Pathogenic",
"Conflicting interpretations of pathogenicity")

gnomad = read_csv("../test_data/Inputs/1-2-3-5_DNMT3A1_gnomad_file.csv",col_names = T) |>
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
         #Gene_name = "DNMT3A1",
         #Prot_ID = "Q9Y6K1") 
  select(Position,Original_Res,New_Res)

arch = read_tsv("../test_data/Inputs/4_DNMT3A1_architecture_file.txt")
ptm = read_tsv("../test_data/Inputs/4_DNMT3A1_post_translation_file.txt") %>% 
  mutate(Annotation="PTM")
plength = 912    

mut = tibble(site=as.numeric(gnomad$Position),Annotation="MissenseVariant")
annotation = bind_rows(mut,ptm) |>
  group_by(Annotation) %>%
  mutate(y = cur_group_id())

arch |> ggplot() +
  geom_rect(xmin = 0,xmax=plength,ymin=-0.2,ymax=0.2,fill="darkgrey")+
  geom_rect(aes(xmin=start_site,xmax=end_site,ymin=-0.5,ymax=0.5),fill=domain_fill)+
  geom_text(aes(x=start_site+((end_site-start_site)/2), y=0, label=architecture_name),size=4,colour="white") +
  theme_bw()+
  scale_x_continuous(breaks = seq(0,plength,25),labels = seq(0,plength,25))+
  theme(axis.text.x = element_text(vjust = 1,angle = 90),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        ) +
  geom_point(data = annotation,mapping=aes(x=site,colour=Annotation,y=y+1),alpha=0.5,size=3)+
  scale_colour_manual(values = cols)+
  coord_cartesian(ylim=c(-2,max(annotation$y)+2),xlim=c(0,plength))+
  labs(x="Position",title = pname)
  
  


