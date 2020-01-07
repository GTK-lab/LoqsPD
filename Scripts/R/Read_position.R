library(here)
library(tidyverse)
library(cowplot)
library(yaml)
library(reshape2)


setwd("/mnt/gtklab01/jinwee/legit_stuff")
# 'Okamura','GSE37443Head','GSE37443Ovary'
dataset <- 'Okamura'
runID <- '21filter_multi'
strand_dir <- '21-21strand_output' #set / hard coded here 
feature <- 'encodecisnat' #encodecisnat
main_dir <- getwd()


#Setting directories
counts_dir <- paste(main_dir,'Output', dataset , 'Mapping', paste('runID=',runID,sep = ''), strand_dir, feature, sep="/")
general_counts_dir <-  paste(main_dir,'Output', dataset , 'Mapping', paste('runID=',runID,sep = ''), 'Counts_output', sep="/")
fig_dir <- paste(main_dir,'Output', 'Figures', sep="/")
config_data <- yaml.load_file(paste(main_dir,'Configfiles','main_config.yaml', sep = '/'))

lib_dict <- config_data[[paste(dataset,'_lib_dict',sep='')]]
lib_list <- names(lib_dict)
seq_index_list <- config_data[[paste(runID, '_index_list', sep ='')]]
non_seq_index_list <- config_data$nonseq_index_list
full_index_list <- append(seq_index_list, non_seq_index_list)
index_dict <- config_data[['index_dict']]

normalization_df <- data.frame(Libraries = names(lib_dict))

for (index in non_seq_index_list){
  norm_file <- file.path(general_counts_dir, paste(index, '_rawcounts.txt',sep=''))
  if(file.exists(norm_file)){
    norm_data <- read_tsv(norm_file)
    normalization_df <- merge(normalization_df, norm_data, by='Libraries')
  }
}

lib_order <- normalization_df$Libraries
rownames(normalization_df) <- normalization_df$Libraries

if (dataset == 'Okamura'){
  norm_type <- "spikeinset01"
  rows <- 2
  dist_color_palette <- c('gray50','goldenrod2','skyblue4','indianred3')
  lab <- c(expression(paste(italic("loqs")^"KO")),"LoqsPA",
           "LoqsPB",
           "LoqsPD")
}else{
  norm_type <- "dmel-all-chromosome-r6.22"
  rows <- 4
  dist_color_palette <- c('moccasin','indianred4','lightcoral','lightpink2',
                          'darkseagreen','goldenrod2','steelblue3','slategray3','slategray','gray45')
  lab <-c(expression(paste(italic("wild-type"))),
          expression(paste(italic("loqs")^"KO ", "; Loqs-PA,PB,PD")),
          expression(paste(italic("loqs")^"KO ", "; Loqs-PA,PD")),
          expression(paste(italic("loqs")^"KO ", "; Loqs-PB,PD")),
          expression(paste(italic("loqs")^"KO ", "; Loqs-PA,PB")),
          expression(paste(italic("loqs")^"KO ", "; Loqs-PA")),
          expression(paste(italic("loqs")^"KO ", "; Loqs-PB")))
}

plus_count_df <- read_tsv(file.path(counts_dir, paste(lib_list[1], feature, 'plus','readpos.txt',sep='_')))
minus_count_df <- read_tsv(file.path(counts_dir, paste(lib_list[1], feature, 'minus','readpos.txt',sep='_')))
plus_count_df$group_counts <- (plus_count_df$group_counts * 1000) / normalization_df[lib_list[1],norm_type]
minus_count_df$group_counts <- (minus_count_df$group_counts * 1000) / normalization_df[lib_list[1],norm_type]

for (lib in lib_list[-1]){
  plus_counts <- read_tsv(file.path(counts_dir, paste(lib, feature, 'plus','readpos.txt',sep='_')))
  minus_counts <- read_tsv(file.path(counts_dir, paste(lib, feature, 'minus','readpos.txt',sep='_')))
  
  plus_counts$group_counts <- (plus_counts$group_counts * 1000) / normalization_df[lib,norm_type]
  minus_counts$group_counts <- (minus_counts$group_counts * 1000) / normalization_df[lib,norm_type]
  
  plus_count_df <- rbind(plus_count_df, plus_counts)
  minus_count_df <- rbind(minus_count_df, minus_counts)
}


TE_feat <- '1731'
if (feature == 'transposable_element'){
  plus_count_df <- plus_count_df %>%
    filter(str_detect(loci, TE_feat))
  minus_count_df <- minus_count_df %>%
    filter(str_detect(loci, TE_feat))
}


##converting $loci to readable
loci_readable <- function(locis, feature){
  if (feature == "transposable_element"){
    output <- sapply(as.character(locis), function(x) strsplit(x, ";")[[1]][2]) %>%
      sapply(function(x) strsplit(x, "=")[[1]][2]) %>%
      unname()
  }else if (feature == "hpRNA_merge"){
    output <- sapply(as.character(locis), function(x) strsplit(x, ";")[[1]][2]) %>%
    sapply(function(x) strsplit(x, "=")[[1]][2]) %>%
    sapply(function (x) strsplit(x, "}")[[1]][1]) %>% 
    unname()
  }else{
    output <- locis
  } 
  output
}

minus_count_df$condition <- sapply(sapply(minus_count_df$lib, function(x) lib_dict[[x]]), function(y) AssignLibCondition(y))
plus_count_df$condition <- sapply(sapply(plus_count_df$lib, function(x) lib_dict[[x]]), function(y) AssignLibCondition(y))

plus_count_df$prime5 <- plus_count_df$sense_pos
plus_count_df$prime3 <- plus_count_df$sense_pos + plus_count_df$length

minus_count_df$prime3 <- minus_count_df$sense_pos
minus_count_df$prime5 <- minus_count_df$sense_pos + minus_count_df$length

plus_count_df$loci <- loci_readable(plus_count_df$loci, feature)
minus_count_df$loci <- loci_readable(minus_count_df$loci, feature)

minus_5prime <- minus_count_df %>%
  group_by(prime5, loci, lib,condition) %>%
  summarise(summed_counts = sum(group_counts))

plus_5prime <- plus_count_df %>%
  group_by(prime5, loci, lib,condition) %>%
  summarise(summed_counts = sum(group_counts))

## Extra filtering
#value <- ifelse(dataset=='Okamura', 0.1, 0.00001) 
#minus_5prime <- minus_5prime[minus_5prime$summed_counts >= value ,]
#plus_5prime <- plus_5prime[plus_5prime$summed_counts >= value,]

if (dataset=='Okamura'){
  mean_minus_5prime <- group_by(minus_5prime, prime5, loci, condition) %>% 
    summarise(mean_counts = mean(summed_counts), sd= sd(summed_counts), se = sd(summed_counts)/sqrt(n()))
  mean_plus_5prime <- group_by(plus_5prime, prime5, loci, condition) %>% 
    summarise(mean_counts = mean(summed_counts), sd= sd(summed_counts), se = sd(summed_counts)/sqrt(n()))
}else{
  mean_minus_5prime <- minus_5prime
  mean_plus_5prime <- plus_5prime
  mean_minus_5prime$mean_counts <- mean_minus_5prime$summed_counts
  mean_plus_5prime$mean_counts <- mean_plus_5prime$summed_counts
}


if (dataset != "Okamura"){
  mean_minus_5prime <- mean_minus_5prime[!mean_minus_5prime$condition %in% c("loqsf00791", "loqsf00791/CyO", "LoqsKO/CyO"),]
}

if (feature == 'hpRNA_merge'){
  test_minus_ft <- c("hpRNA:CR46342", "hpRNA:CR18854")
} else if (feature == 'encodecisnat'){
  test_minus_ft <- c("chr3L:25114918-25115185", "chr3L:15560679-15561111")
}

test_mean_minus5 <- mean_minus_5prime[mean_minus_5prime$loci %in% test_minus_ft, ] 

mean_minus5_plot <- ggplot(test_mean_minus5, aes(x= prime5,y= mean_counts,color=condition)) +
  geom_step(size=0.2)+
  scale_colour_manual(values=dist_color_palette,
                      labels=lab) +
  #facet_grid(condition ~ loci, scales = "free") +
  facet_grid(loci ~ condition, scales = "free") +
  labs(colour="Rescue Condition",
       y="Mean Counts per thousand spike-in",
       x = "5' Start Position") +
  theme(axis.text.x=element_text(size = 8),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10),
        axis.text.y=element_text(size = 8),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        #strip.background = element_blank(),
        panel.background = element_blank(),
        #$rect = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        #
        #strip.text.x = element_blank(),
        # legend.title = element_text(size=6),
        # legend.text = element_text(size=4),
        legend.position = "none"
  )

mean_minus5_plot

save_plot(file.path(fig_dir,paste(dataset, feature,'minusreadpos.pdf',sep='-')),
          mean_minus5_plot,
          base_aspect_ratio = 2)