#Input files: {lib}_counts.txt, main_config.yaml
library(here) #here() starts at /mnt/raid0/home/jinwee
library(tidyverse)
library(cowplot)
library(yaml)
library(gtools)
#library(devtools)
library(superheat)
library(reshape)
library(reshape2)
library(cetcolor)
library(xlsx)
#library(edgeR)
##'Okamura','GSE37443Head','GSE37443Ovary'
#test variables
setwd("/mnt/gtklab01/jinwee/legit_stuff")
dataset <- 'GSE37443Ovary'
runID <- '21filter'
FC_feature <-'transposable_element'
main_dir <- getwd()

counts_dir <- paste(main_dir,'Output', dataset , 'Mapping', paste('runID=',runID,sep = ''), 'Counts_output', sep="/")
fc_counts_dir <- paste(main_dir,'Output', dataset , 'Mapping', paste('runID=',runID,sep = ''), 'FC_output', paste(FC_feature,'FC',sep='_'), sep="/")
fig_dir <- paste(main_dir, 'Output', 'Figures', sep='/')

config_data <- yaml.load_file(paste(main_dir,'Configfiles','main_config.yaml', sep = '/'))

#Load config data, lib_dict and merge both seq and non_seq indexes 
lib_dict <- config_data[[paste(dataset,'_lib_dict',sep='')]]
lib_list <- names(lib_dict)
non_seq_index_list <- config_data$nonseq_index_list

count_df <- read_tsv(file.path(fc_counts_dir, paste(lib_list[1], FC_feature, 'locicounts.txt',sep='_')))

for (lib in lib_list[-1]){
  lib_counts <- read_tsv(file.path(fc_counts_dir, paste(lib, FC_feature, 'locicounts.txt',sep='_')))
  count_df <- smartbind(count_df,lib_counts,fill=0)
}

#Casting entire DE_df to numeric because smartbind converts to character 
count_df <- as.data.frame(apply(count_df, 2, as.numeric, simplify=TRUE))
count_df <- count_df + 1 

#Adding rownames and transposing
rownames(count_df) <- lib_list

count_df <- t(count_df) %>%
  as.data.frame()

check_df <- count_df
feature_order <- rownames(count_df)

#Setting up for normalization and implementing normalization 
normalization_df <- data.frame(Libraries = names(lib_dict))

for (index in non_seq_index_list){
  count_data <- read_tsv(file.path(counts_dir, paste(index, '_rawcounts.txt',sep='')))
  normalization_df <- merge(normalization_df, count_data, by='Libraries')
}

lib_order <- normalization_df$Libraries

if(dataset=='Okamura'){
  normalization_df <- subset(normalization_df, select= c('spikeinset01'))
  rownames(normalization_df) <- lib_order
  norm_count_df <- sapply(lib_list, function (x) sapply(as.numeric(as.vector(count_df[,x])), function(y) y * 1000 /normalization_df[x,])) %>%
    as.data.frame()
}else{
  normalization_df <- subset(normalization_df, select= c('dmel-all-chromosome-r6.22'))
  rownames(normalization_df) <- lib_order
  norm_count_df <- sapply(lib_list, function (x) sapply(as.numeric(as.vector(count_df[,x])), function(y) y * 1000000 /normalization_df[x,])) %>%
    as.data.frame()
}

#Ok now just adding simple barplot 
rownames(norm_count_df) <- feature_order

AssignLibCondition <- function(x) {
  if (grepl('Control', x) == TRUE){
    'Control'
  }else if (grepl('_rep', x)){
    strsplit(x, '_rep')[[1]][1]
  }else{
    x
  }
}


norm_count_df <- t(norm_count_df) %>%
  as.data.frame()

Condition <- as.vector(sapply(rownames(norm_count_df), FUN = AssignLibCondition))

norm_count_df$Condition <- Condition
norm_count_df$Libraries <- rownames(norm_count_df)


melted_normcount_df <- subset(norm_count_df, select = -c(Condition)) %>%
  melt(id.vars = 'Libraries',
       measure.vars = feature_order,
       variable.name = 'feature',
       value.name = 'normcounts')

if (FC_feature == 'hpRNA_merge' | FC_feature == 'transposable_element'){
  melted_normcount_df$feature <- sapply(as.character(melted_normcount_df$feature), function(x) strsplit(x, ";")[[1]][2]) %>%
    sapply(function(x) strsplit(x, "=")[[1]][2]) %>%
    sapply(function (x) strsplit(x, "}")[[1]][1])
  melted_normcount_df$feature <- substr(as.character(melted_normcount_df$feature), 1,nchar(as.character(melted_normcount_df$feature))-1)
}


melted_sumcount_df <- melted_normcount_df %>% 
  group_by(Libraries, feature) %>%
  summarise(sum_counts = sum(normcounts))

melted_sumcount_df$Libraries <- sapply(melted_sumcount_df$Libraries, function (x) lib_dict[[x]])
melted_sumcount_df$Condition <- sapply(melted_sumcount_df$Libraries, FUN = AssignLibCondition)
melted_sumcount_df <- melted_sumcount_df[order(melted_sumcount_df$Condition),]

melted_mean_df <- group_by(melted_sumcount_df,Condition,feature) %>%
  summarise(mean_sumcounts=mean(sum_counts), sd=sd(sum_counts), se = sd(sum_counts)/sqrt(n()))


if (dataset=='Okamura'){
  counts_y_label <- 'Counts per thousand spike-in (+/- standard error)'
  counts_x_label <- 'Rescue isoform'
  color_palette <- c('gray45', 'goldenrod2','skyblue4','indianred3')
  ytitle_size <- 15
  xtitle_size <- 13
  xtext_size <- 4.5
  ytext_size <- 8.5
  patch_strip_size <- 6
  patch_xtext_size <- 5
  patch_label_size <- 10 ## default 10
  mylabels <- c(expression(paste(italic("loqs")^"KO")),
                "LoqsPA",
                "LoqsPB",
                "LoqsPD")
}else{ 
  color_palette <- c('moccasin',
                     'indianred4',
                     'lightcoral',
                     'lightpink2',
                     'darkseagreen',
                     'goldenrod2',
                     'steelblue3',
                     'slategray3',
                     'slategray',
                     'gray45'
  )
  patch_xtext_size <- 4
  patch_ytext_size <- 5
  patch_strip_size <- 6
  patch_label_size <- 12
  patch_title_size <- 2
  patch_vjust <- 0.5
  counts_y_label <- 'Counts per million'
  counts_x_label <- 'Library condition'
  library_order <- c('Control',
                     'LoqsPA+LoqsPB+LoqsPD',
                     'LoqsPA+LoqsPD',
                     'LoqsPB+LoqsPD',
                     'LoqsPA+LoqsPB',
                     'LoqsPA',
                     'LoqsPB',
                     'loqsf00791/CyO',
                     'loqsf00791',
                     'LoqsKO/CyO'
  )
  mylabels <- c(expression(paste(italic("wild-type"))),
                expression(paste(italic("loqs")^"KO ", "; Loqs-PA,PB,PD")),
                expression(paste(italic("loqs")^"KO ", "; Loqs-PA,PD")),
                expression(paste(italic("loqs")^"KO ", "; Loqs-PB,PD")),
                expression(paste(italic("loqs")^"KO ", "; Loqs-PA,PB")),
                expression(paste(italic("loqs")^"KO ", "; Loqs-PA")),
                expression(paste(italic("loqs")^"KO ", "; Loqs-PB")),
                expression(paste(italic("f00791"), "/CyO")),
                expression(paste(italic("f00791"))),
                expression(paste(italic("loqs")^"KO", "/CyO"))
  )
  melted_mean_df$Condition <- factor(melted_mean_df$Condition, levels=library_order)
  melted_mean_df <- melted_mean_df[order(melted_mean_df$Condition),]
}

  
if (dataset == 'Okamura'){
  prior_lst1 <-  c('mdg1', '1731', '297', 'roo', 'springer')
  prior_lst2 <-  c('gypsy3', 'gypsy6', 'Doc', 'blood', 'S')
}else if (dataset == 'GSE37443Head'){
  prior_lst1 <-  c('1731', 'gypsy3', 'gypsy6','blood', 'S')
  prior_lst2 <-  c('mdg1','Doc', '297', 'roo', 'springer')
}else{
  prior_lst1 <-  c('gypsy3', 'gypsy6', 'blood','mdg1', 'Stalker2', 'Stalker4' )
  prior_lst2 <-  c('roo','springer','F', 'Doc', '297','412')
}

siRNA_mean_df1 <- melted_mean_df[melted_mean_df$feature %in% prior_lst1,]
siRNA_mean_df2 <- melted_mean_df[melted_mean_df$feature %in% prior_lst2,]
siRNA_mean_df1$foo <- 1
siRNA_mean_df2$foo <- 2

full_siRNA_df <- rbind(siRNA_mean_df1,siRNA_mean_df2)
  
TE_siRNA_plot <- ggplot(full_siRNA_df, aes(x= feature, y=mean_sumcounts, fill =Condition)) +
    geom_bar(stat="identity",colour="black", size=0.3 , position=position_dodge())+
    facet_wrap( ~ foo, nrow = 2, ncol=1, scales = "free") +
    geom_errorbar(aes(ymin=mean_sumcounts-se, ymax=mean_sumcounts+se),
    size= 0.3,
    width=.2,
    position=position_dodge(0.9)) + 
    labs(y = counts_y_label, x= "TE Feature") +
    scale_fill_manual(values=color_palette,
                      labels= mylabels,
                      name = "Rescue Isoform")+
  theme_bw()+
  theme(
      axis.title.x = element_text(size=9),
      axis.text.x = element_text(size=7),
      axis.title.y = element_text(size=9),
      axis.text.y = element_text(size=8),
      legend.text = element_text(size = 4),
      legend.title = element_text(size = 6),
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      rect = element_blank(),
      axis.line.x = element_line(color="black"),
      axis.line.y = element_line(color="black")
    ) 

TE_siRNA_plot
  
save_plot(file.path(fig_dir,paste(dataset,'TEcounts.pdf',sep='-')),
            TE_siRNA_plot)
