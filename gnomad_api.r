
querymaker <- function(gene,dataset_id){
  paste0('{\n  gene(gene_symbol: "', gene, '", reference_genome: GRCh37) {
    variants(dataset: ',dataset_id,') {
      variant_id
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
}')}


query <- querymaker(gene = "MECP2", dataset_id = "gnomad_r2_1")

library(httr2)
library(tidyverse)
jsondata <- request("https://gnomad.broadinstitute.org/api/?") %>%
  req_body_json(list(query=query, variables="null")) %>%
  req_perform() %>%
  resp_body_json()

gnomad_v<-jsondata$data$gene$variants |> map(~unlist(.x) |> t() |> as.data.frame()) |> bind_rows()
gnomad_c<-jsondata$data$gene$clinvar_variants |> map(~unlist(.x) |> t() |> as.data.frame()) |> bind_rows()

gnomad <- gnomad_v |> left_join(gnomad_c,by="variant_id")

gnomad_v |> left_join(gnomad_c,by="variant_id") |>
  mutate(across(starts_with(c("genome","exome")),~coalesce(as.numeric(.x),0))) |>
  mutate(af = (genome.ac+exome.ac) / (genome.an+exome.an)) |>
  select(`VEP Annotation` = consequence,
         `Protein Consequence` = hgvsp,
         `ClinVar Clinical Significance` = clinical_significance,
         `Allele Frequency` = af)

dim(gnomad)
dim(gnomad_v)
