library(devtools)
## devtools::install_github("thomasp85/patchwork")
library(ggplot2)
library(patchwork)
library(magick)

setwd("/mnt/gtklab01/jinwee/legit_stuff")
main_dir <- getwd()
fig_dir <- paste(main_dir,'Output', 'Figures', sep="/") ## note create this directory first

pdf_flowchart <- "/mnt/raid0/home/jinwee/loqs_proj_v2/Output/patchwork/flowchart.pdf" 


##Fig1
flowchartpdf <- image_read_pdf(pdf_flowchart)
flowchart <- image_ggplot(flowchartpdf)
flowchart <- flowchart + labs(tag = "A")

okamura <- patch_normcount_plot + patch_FC_plot + plot_layout(ncol = 1)
fig1 <- flowchart + okamura
save_plot(file.path(dir,'Figure-1.pdf'),
          fig1,
          base_aspect_ratio = 2)


##Fig2 
fc_plots <- patch_FC_plot  + patch_FC_plot2
fig2 <- patch_normcount_plot + patch_normcount_plot2 - fc_plots + plot_layout(ncol = 1)

save_plot(file.path(fig_dir,"Figure-2.pdf"),
          fig2,
          base_aspect_ratio = 2)
  






