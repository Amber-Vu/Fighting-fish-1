#  GoM
source("https://bioconductor.org/biocLite.R")
biocLite("cellTree")

library(devtools)
install_github("kkdey/maptpx")
install_github('kkdey/CountClust')

library(maptpx)
library(cellTree)
library(CountClust)


data <- read.table("input file", header=TRUE, row.names=1, quote="")

lda.results <- compute.lda(data,k.topics=7:7, method="maptpx")

annotation <- data.frame(
  sample_id = paste0("X",c(1:NROW(lda.results$omega))),
  tissue_label = factor(rownames(lda.results$omega),
                        levels = rev(c("B_1","B_2","B_3","B_4","B_5",
                                       "D2_11","D2_12","D2_21","D2_22","D2_31", "D2_32", "D2_41", "D2_42", "D2_51", "D2_52",
                                       "D6_11","D6_12", "D6_21", "D6_22", "D6_31","D6_32", "D6_41", "D6_42", "D6_51", "D6_52"))))


annotation
StructureGGplot(omega=lda.results$omega, annotation=annotation, palette = RColorBrewer::brewer.pal(8, "Accent"),
                order_sample=TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))



