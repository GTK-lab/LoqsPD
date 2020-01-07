#Input files: {lib}_counts.txt, main_config.yaml
library(here) #here() starts at /mnt/raid0/home/jinwee
library(tidyverse)
library(cowplot)
library(yaml)
library(gtools)
#library(devtools)
library(superheat)
library(reshape2)
library(cetcolor)
#library(edgeR)
## Okamura','GSE37443Head','GSE37443Ovary', 'lengthfilter'
##'lengthfilter'
#test variables 

setwd("/mnt/gtklab01/jinwee/legit_stuff")

dataset <- 'GSE37443Ovary'
runID <- '21filter'
FC_feature <-'hpRNA_merge'
main_dir <- getwd()
#['antiexon',03_kawamura','01_okamura', '02_czech', '04_ghildiyal','hpRNA_merge','transposable_element']

counts_dir <- paste(main_dir,'Output', dataset , 'Mapping', paste('runID=',runID,sep = ''), 'Counts_output', sep="/")
fc_counts_dir <-  paste(main_dir,'Output', dataset , 'Mapping', paste('runID=',runID,sep = ''), 'FC_output', paste(FC_feature, 'FC', sep='_'), sep="/")
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

rownames(norm_count_df) <- sapply(rownames(norm_count_df), function (x) lib_dict[[x]])
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
    sapply(function(x) strsplit(x, "=")[[1]][2])
}


melted_normcount_df$Condition <- sapply(melted_normcount_df$Libraries, FUN = AssignLibCondition)
melted_normcount_df <- melted_normcount_df[order(melted_normcount_df$Condition),]


if (FC_feature == 'antiexon'){
  melted_normcount_df$variable <- process_antiexon(as.character(melted_normcount_df$variable))
  melted_normcount_df <- melted_normcount_df %>%
    group_by(Libraries, variable) %>% 
    summarise(value = sum(value))
  melted_normcount_df$Condition <- sapply(melted_normcount_df$Libraries, FUN = AssignLibCondition)
}


melted_FC_normcount_df <- group_by(melted_normcount_df,Condition,feature) %>%
  summarise(mean_normcounts=mean(normcounts), sd=sd(normcounts), se = sd(normcounts)/sqrt(n()))


# 
# ##saving
# save_mean_df <-dcast(melted_FC_normcount_df, variable~Condition, value.var = "mean_normcounts")
# rownames(save_mean_df) <- save_mean_df$variable
# save_mean_df <- select(save_mean_df, -c("variable"))
# save_mean_df <- save_mean_df[rev(order(rowSums(save_mean_df))),]
# write.xlsx(save_mean_df, file= here("loqs_proj_v2","okamura_exonantisense.xlsx"), sheetName="mean_norm_counts", row.names=TRUE, col.names = TRUE, append = TRUE)
# 


if (dataset == "Okamura"){
  melted_FC_normcount_df$feature <- factor(melted_FC_normcount_df$feature,
                                           levels=c("hpRNA:CR46342","hpRNA:CR18854","hpRNA:CR32205",
                                                    "hpRNA:CR32207","hpRNA:CR33940", "hpRNA:1", "hpRNA:CR46340"))
} else if (dataset == "GSE37443Head"){
  melted_FC_normcount_df$feature <- factor(melted_FC_normcount_df$feature,
                                           levels=c("hpRNA:CR46342","hpRNA:CR46343","hpRNA:CR18854","hpRNA:CR32205",
                                                    "hpRNA:CR32207","hpRNA:CR33940", "hpRNA:1"))
  melted_FC_normcount_df <- filter(melted_FC_normcount_df, feature != 'hpRNA:1') %>% 
    filter(feature != 'hpRNA:CR46343') %>% 
    filter(feature != 'hpRNA:CR32205')
} else if (dataset == "GSE37443Ovary"){
  melted_FC_normcount_df$feature <- factor(melted_FC_normcount_df$feature,
                                           levels=c("hpRNA:CR46342","hpRNA:CR18854","hpRNA:CR32205",
                                                    "hpRNA:CR32207","hpRNA:CR33940", "hpRNA:1"))
  melted_FC_normcount_df <- filter(melted_FC_normcount_df, feature != 'hpRNA:1') %>% 
    filter(feature != 'hpRNA:CR32207') %>%
    filter(feature != 'hpRNA:CR32205')
}

melted_FC_normcount_df <- melted_FC_normcount_df[order(melted_FC_normcount_df$feature),] 



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
  final_features <- c("hpRNA:CR46342","hpRNA:CR18854","hpRNA:CR32205","hpRNA:CR32207","hpRNA:CR33940")
  melted_FC_normcount_df <- melted_FC_normcount_df[melted_FC_normcount_df$variable %in% final_features,]
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
  library_order <- c('WT_Control',
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
  melted_FC_normcount_df[melted_FC_normcount_df=='Control'] <- 'WT_Control'
  melted_FC_normcount_df$Condition <- factor(melted_FC_normcount_df$Condition, levels=library_order)
  melted_FC_normcount_df <- melted_FC_normcount_df[order(melted_FC_normcount_df$Condition),]
  final_features <- unique(melted_FC_normcount_df$feature)
}

newline_features <- gsub(":", ":\n", final_features)
names(newline_features) <- final_features

patch_FC_plot2 <- ggplot(melted_FC_normcount_df, aes(x= Condition, y=mean_normcounts)) +
  geom_bar(stat="identity",colour="black", size=0.3 ,aes(fill =Condition))+
  scale_x_discrete(labels = mylabels) +
  # geom_errorbar(aes(ymin=mean_normcounts-se, ymax=mean_normcounts+se),
  # size= 0.3,
  # width=.2,
  # position=position_dodge(.9))+
  facet_wrap( ~ feature, nrow=1, scales = "free_x", labeller= as_labeller(newline_features)) +
  coord_flip() +
  scale_fill_manual(values=color_palette) +
  #scale_fill_brewer(palette="Spectral") +
  labs(y = counts_y_label, x=counts_x_label) + 
  #labs(title = paste('dataset',dataset,sep='='),subtitle = paste('runID',runID,sep='='))+
  #theme(plot.title = element_text(hjust=0, size = 9),plot.subtitle = element_text(size = 8))+
  #theme(strip.text.x =element_text(size = 8, margin = margin(0.1, 0, 0.1, 0, "cm"))) +
  theme(legend.position ="None",legend.text = element_text(size = 5), legend.title = element_text(size = 6)) + 
  theme(axis.text.x=element_text(size = patch_xtext_size,angle=30)) +
  theme(axis.text.y=element_text(size = patch_ytext_size)) +
  theme(axis.title=element_text(size=patch_label_size)) + 
  ##labs(tag = "C") +
  ##theme(axis.title.x = element_text(size =xtitle_size), axis.title.y = element_text(size =ytitle_size)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        rect = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        ##axis.title.x=element_blank(),
        ##axis.title.y=element_blank(),
        strip.text.x = element_text(size=patch_strip_size, face="bold", hjust = 0.1, vjust=patch_vjust))


patch_FC_plot2

