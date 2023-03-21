#This code is for plotting VdVp plots from data downloaded from the VdVp Plotting Tool

#fetches library
library("ggplot2")

#defines parameters set by user
protein_name <- "DNMT3B"
protein_length <- 853
x_axis_breaks <- 50
path_to_data <- "/Users/gaurideak/Desktop/testing/Plot_Files/input.txt"
path_for_export <- "/Users/gaurideak/Desktop/testing/output.tiff"

#loads input file downloaded from VdVp Plotting Tool in txt format
df <- read.delim(file=path_to_data, sep=",",
                 header=TRUE)
#smooths data
splined_data <- as.data.frame(spline(df$Residue, df$VdVp))

#plots data as line graph
p <- ggplot() +
  theme_classic() +
  geom_line(data = splined_data, aes(x = x, y = y), size = 0.4) +
  scale_x_continuous(breaks=seq(0,protein_length,x_axis_breaks)) +
  labs(title = protein_name, x = "Residue", y = "Vd/Vp") +
  theme(plot.title = element_text(size=16, face="bold"),
        axis.text=element_text(size=7, color = "black"),
        axis.title=element_text(size=11, face="bold"))
p

#save plot with custom dimensions
ggsave(file=path_for_export, width=6, height=2, dpi=600)

