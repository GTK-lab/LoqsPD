#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

library(here)
library(tidyverse)
library(cowplot)
library(yaml)
library(reshape2)
library(dplyr)

#external variables for Rscript
#dataset <- args[1]
#runID <- args[2]

setwd("/mnt/gtklab01/jinwee/legit_stuff")

# 'Okamura','GSE37443Head','GSE37443Ovary'
dataset <- 'GSE37443Head'
runID <- '21filter'
main_dir <- getwd()
#Setting directories
##counts_dir <- here('loqs_proj_v2','Output', dataset , 'Mapping', paste('runID=',runID,sep = ''),'Counts_output')
counts_dir <- paste(main_dir,'Output', dataset , 'Mapping', paste('runID=',runID,sep = ''), 'Counts_output', sep="/")
# figures_dir <- here('loqs_proj_v2','Output', dataset , 'Mapping', paste('runID=',runID,sep = ''),'Figures_output')
# maplog_dir <- here('loqs_proj_v2','Output', dataset , 'Mapping', 'Maplog')

#Load config data, lib_dict and merge both seq and non_seq indexes 
config_data <- yaml.load_file(paste(main_dir,'Configfiles','main_config.yaml', sep = '/'))

lib_dict <- config_data[[paste(dataset,'_lib_dict',sep='')]]

seq_index_list <- config_data[[paste(runID, '_index_list', sep ='')]]
non_seq_index_list <- config_data$nonseq_index_list
full_index_list <- append(seq_index_list, non_seq_index_list)
full_index_list <- setdiff(full_index_list,
                           c("rRNARM","rRNA_exon_merge","chr1dmerdna","RNARM","tRNA_exon_merge","snRNA_exon_merge","snoRNA_exon_merge"))
full_index_list
index_dict <- config_data[['index_dict']]


#Initialize base count_df here 
#note change the library tags to type later
count_df <- data.frame(Libraries = names(lib_dict))

#Load data into df row by row 
for (index in full_index_list){
  count_data <- read_tsv(file.path(counts_dir, paste(index, '_rawcounts.txt',sep='')))
  count_df <- merge(count_df, count_data, by='Libraries')
}

#Normalize by RPM and spikein if dataset='Okamura', store Library column as variable to prevent loss
#when performing numerical operation on entire dataframe. 
lib_column <- count_df$Libraries

count_df <- subset(count_df, select =-c(Libraries)) + 1 
count_df$Libraries <- lib_column

if(dataset=='Okamura'){
  norm_count_df <- subset(count_df, select = -c(spikeinset01,`dmel-all-chromosome-r6.22`,Libraries)) * 1000 / count_df$spikeinset01
}else{
  norm_count_df <- subset(count_df, select = -c(spikeinset01,`dmel-all-chromosome-r6.22`,Libraries)) * 1000000 / count_df$`dmel-all-chromosome-r6.22`
}

#Collapse + sum columns/indexes based on their general index_type
colnames(norm_count_df) <- sapply(colnames(norm_count_df), function(x) index_dict[[x]])
all_index_types <- colnames(norm_count_df)
uniq_index_types <- unique(colnames(norm_count_df)) #Providing some sort of reference to deal with repeated column names 

#Initialize norm_count_df_col here
norm_count_df_col <- data.frame(data.frame(row.names = 1:nrow(norm_count_df))) 

for (index_type in uniq_index_types){
  index_type_cols <-as_data_frame(norm_count_df[, all_index_types ==index_type])
  norm_count_df_col[index_type] <- rowSums(index_type_cols)
}


#Ok, lets try removing Background, ncRNA, pseudo, remove ncRNA and Pseudogene from the diagram
#and make it a 1 row thing 
#norm_count_df_col <- select(norm_count_df_col, -c(Background))
colnames(norm_count_df_col) <- c("miRNA", "hpRNA","TE-derived", "cis-NAT", "other", "exonic-antisense")
final_index_types <- colnames(norm_count_df_col)
  
#Function to label Library_type with Condition, '_rep' marks a replicate library, standardized naming conventions
AssignLibCondition <- function(x) {
  if (grepl('Control', x) == TRUE){
    'Control'
  }else if (grepl('_rep', x)){
    strsplit(x, '_rep')[[1]][1]
  }else{
    x
  }
}


#Reassign values in Libraries column as 'Library_type'
norm_count_df_col$Libraries <- sapply(lib_column, function(x) lib_dict[[x]])
norm_count_df_col$Condition <- sapply(norm_count_df_col$Libraries, FUN = AssignLibCondition)

melted_normcount_df <- subset(norm_count_df_col, select = -c(Condition)) %>%
  melt(id.vars = 'Libraries',
       measure.vars = final_index_types,
       variable.name = 'Index',
       value.name = 'value')

melted_normcount_df$Condition <- sapply(melted_normcount_df$Libraries, FUN = AssignLibCondition)
melted_normcount_df <- melted_normcount_df[order(melted_normcount_df$Condition),]

test <- melted_normcount_df

###Starting updated normcounts plots
col_order <- c('miRNA','hpRNA', 'TE-derived','cis-NAT', 'other','exonic-antisense')
test$variable <- factor(test$Index, levels=col_order)
melted_normcount_df <- test[order(test$variable),]
#c('gray50', 'lightgoldenrod1','skyblue3','lightcoral')

#Figure variables dependant on dataset
if (dataset=='Okamura'){
  counts_y_label <- 'Counts per thousand spike-in (+/- standard error)'
  counts_x_label <- 'Rescue isoform'
  title <- "                    endogenous siRNA                       "
  color_palette <- c('gray45', 'goldenrod2','skyblue4','indianred3')
  ytitle_size <- 17
  xtitle_size <- 1
  xtext_size <- 5.5
  patch_xtext_size <- 5
  patch_ytext_size <- 5
  patch_strip_size <- 6 ## 7
  patch_label_size <- 10 ## default 10
  patch_vjust <- 0.5
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
  counts_y_label <- 'Counts per million'
  counts_x_label <- 'Library condition'
  ytitle_size <- 13
  xtitle_size <- 12
  xtext_size <- 6
  patch_xtext_size <- 4
  patch_ytext_size <- 5
  patch_strip_size <- 4
  patch_label_size <- 12
  patch_title_size <- 15
  patch_vjust <- 0.5
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
                
}

if (dataset!='Okamura'){
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
  melted_normcount_df[melted_normcount_df=='Control'] <- 'WT_Control'
  melted_normcount_df$Condition <- factor(melted_normcount_df$Condition, levels=library_order)
  melted_normcount_df <- melted_normcount_df[order(melted_normcount_df$Condition),]
  }


melted_group_normcount_df <- group_by(melted_normcount_df,Condition,variable) %>%
  summarise(mean_normcounts=mean(value), sd=sd(value), se = sd(value)/sqrt(n()))

melted_group_normcount_df[is.na(melted_group_normcount_df)] <- 0


newline_indexes <- gsub("-", "-\n", final_index_types)
names(newline_indexes) <- final_index_types
newline_indexes[["cis-NAT"]] <- "cis-NAT"

patch_normcount_plot <- ggplot(melted_group_normcount_df, aes(x= Condition, y=mean_normcounts)) +
  ggtitle(expression(underline("    endogenous siRNA     "))) + ##for fukunaga 
  ##ggtitle(expression(underline("       endogenous siRNA       "))) + ##for okamura supp 
  ##ggtitle(expression(underline("           endogenous siRNA          "))) + ## for okamura
  #ggtitle(expression(underline("                    endogenous siRNA                       "))) +
  scale_x_discrete(labels = mylabels) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  geom_bar(stat="identity",colour="black", size=0.3 ,aes(fill =Condition))+
  # geom_errorbar(aes(ymin=mean_normcounts-se, ymax=mean_normcounts+se),
  # size= 0.3,
  # width=.2,
  # position=position_dodge(.9))+
  facet_wrap( ~ variable, nrow=1, scales = "free_x", labeller = as_labeller(newline_indexes)) +
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
  theme(plot.title = element_text(size = 20)) +
  theme(axis.title=element_text(size=patch_label_size)) + 
  labs(tag = "A") +
  ##labs(tag = "B") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        rect = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        strip.text.x = element_text(size=patch_strip_size, face="bold", hjust = 0.1, vjust= patch_vjust),
        axis.title.x=element_blank(),
        plot.title = element_text(face="plain", hjust=0.8,size = 13)) ## default is 14 

patch_normcount_plot

## Final compilation of plots happens in patchwork