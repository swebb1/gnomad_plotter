library("httr2")
library("jsonlite")

seq = "MVHLTPEEKSAVTALWGKVNVDEVGG"

myseq = paste0(">myseq\n",seq)


req <- request("http://protein.bio.unipd.it/fellsws/submit") |> 
  req_headers("Content-Type" = "multipart/form-data") |>
  req_body_multipart(sequence = myseq)

req |> req_dry_run()
res <- req_perform(req)

rb = res |> resp_body_json()
job = rb$jobid

req2 <- request(paste0("http://protein.bio.unipd.it/fellsws/status/",job))
res2 = req_perform(req2)
rb2 = res2 |> resp_body_json()
rb2
result = rb2$names[[1]][[2]]

req3 <- request(paste0("http://protein.bio.unipd.it/fellsws/result/",result))
res3 = req_perform(req3)
rb3 = res3 |> resp_body_json()

hca = rb3$hca |> unlist()
dis = rb3$p_dis |> unlist()

phelix = rb3$p_h |> unlist()
pstrand = rb3$p_e |> unlist()
pcoil = rb3$p_c |> unlist()

df <- data.frame(index = seq(1,str_length(seq),1),
                 Helix=as.numeric(phelix),
                 Strand=as.numeric(pstrand),
                 Coil=as.numeric(pcoil)) |>
  mutate(Coil = Coil * -1) |>
  pivot_longer(-index)

hsc_cols = c(Helix="#99004D",Strand="#FF9900",Coil="#878787")

df |> ggplot(aes(index,value,fill=name))+
  geom_col(position = "identity",alpha=0.5) +
  scale_fill_manual(values = hsc_cols)+
  theme_minimal()

