'''
Run as:

snakemake -s Okamuralib_Sequential_mapping.Snakefile --config runID ='' lengthfilter = '' minlength = '' maxlength= '' 


- 'runID': String which corresponds to the Sequential mapping order detailed in `Configfiles/main_config.yaml`, i.e runID_index_list 
- 'lengthfilter': Boolean true/false that determines whether length filtering is performed 
- 'minlength': Integer which sets the minimum (greater than or equal to) length for a mapped read to be recognized as a count. 
- 'minlength': Integer which sets the maximum (less than or equal to) length for a mapped read to be recognized as a count. 

'''

import os 
from os.path import join as pjoin
from Scripts.Python.helpers import MakeGenomeMapLog, SortMapFiles, GatherBowtieCounts, GatherLengthfilterBowtieCounts
import pandas as pd
import sys

configfile: srcdir("Configfiles/main_config.yaml")

#Defining --config variables
dataset = 'Okamura'
runID  = config['runID'] #default now set as edgeR, for 21 filter? hmnmn ok we set '21filter'
lengthfilter = config['lengthfilter']

if lengthfilter == 'true' :
	minlength = config['minlength']
	maxlength = config['maxlength']


#Setting index variables and directories
bt_index_dir = 'Output/Index_generation/Bowtie_indexes'
seq_index_list = config[runID + '_index_list']
nonseq_index_list = config['nonseq_index_list']
full_index_list = seq_index_list + nonseq_index_list
siRNA_index_list = config['siRNA_index_list'] 
index_dict = config['index_dict']
index_type_list = list(set(index_dict.get(x) for x in seq_index_list))
index_type_list.append('Combined')

#Genome version 
genome_vers = config['genome_vers']

#Setting lib_list based on {dataset}
lib_list = list(config[dataset + '_lib_dict'].keys())

#Setting paths, if they have not been created first rule will create them. 
lib_dir = 'Output/' + dataset + '/Preprocessed_libs/Final_preprocessed_libs'
bowtie_out_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/Bowtie_output'
counts_out_dir = 'Output/' + dataset + '/Mapping/runID=' +  runID +'/Counts_output'
figures_out_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/Figures_output'
maplog_out_dir = 'Output/' + dataset + '/Mapping/Maplog'

rule all:
	input:
		#expand(pjoin(bowtie_out_dir,'{lib}_{nonseq_index}_mapped.txt'), lib = lib_list, nonseq_index = nonseq_index_list)
		#pjoin('Flags', dataset + '_runID=' + runID + '_run_dirs_initialized.done'),
		#expand(pjoin(maplog_out_dir, '{lib}_maplog.txt'), lib = lib_list),
		expand(pjoin(bowtie_out_dir, '{lib}_{seq_index}_mapped.txt'), lib=lib_list, seq_index = seq_index_list),
		expand(pjoin(counts_out_dir, '{index}_rawcounts.txt'), index = full_index_list)
		
'''
rule Make_maprun_dirs:
	output:
		touch(pjoin('Flags', dataset + '_runID=' + runID + '_run_dirs_initialized.done'))
	run:
		shell('snakemake -s Initialize_directories.Snakefile --nolock --config dataset={dataset} runID={runID}') #note the nolock here 
'''

rule Spikein_genome_mapping:
	input:
		in_libs = pjoin(lib_dir,'{lib}_final_col_ntfilter.fasta'),
		in_fw_indexes = expand(pjoin(bt_index_dir,'{nonseq_index}.{num}.ebwt'), nonseq_index = nonseq_index_list, num = [1,2,3,4]),
		in_rev_indexes = expand(pjoin(bt_index_dir,'{nonseq_index}.rev.{num}.ebwt'), nonseq_index = nonseq_index_list, num = [1,2])
	output:
		single_mapped = pjoin(bowtie_out_dir,'{lib}_{nonseq_index}_mapped.txt'),
		#multi_mapped = pjoin(bowtie_out_dir,'{lib}_{nonseq_index}-multi_mapped.txt'),
		#single_unmapped = pjoin(bowtie_out_dir,'{lib}_{nonseq_index}_unmapped.fasta'),
		#multi_unmapped = pjoin(bowtie_out_dir,'{lib}_{nonseq_index}-multi_unmapped.fasta')	
	run:
		index_path = pjoin(bt_index_dir, wildcards.nonseq_index)
		#shell('bowtie -f --un {output.multi_unmapped} -a -v 0 -p 1 -t {index_path} {input.in_libs} > {output.multi_mapped}')
		shell('bowtie -f -v 0 -p 1 -t {index_path} {input.in_libs} > {output.single_mapped}')

'''
rule Generate_map_maplog:
	input:
		in_genome_map = pjoin(bowtie_out_dir, '{lib}_' + genome_vers + '-multi_mapped.txt'), 
		in_genome_umapped = pjoin(bowtie_out_dir, '{lib}_' + genome_vers + '-multi_unmapped.fasta')
	output:
		out_maplog = pjoin(maplog_out_dir, '{lib}_maplog.txt')
	run:
		 MakeGenomeMapLog(input.in_genome_map,
		 				  input.in_genome_umapped,
		 				  wildcards.lib,
		 				  output.out_maplog)
'''

rule Sequential_mapping:
	input:
		in_libs = expand(pjoin(lib_dir,'{lib}_final_col_ntfilter.fasta'), lib=lib_list),
		in_fw_indexes = expand(pjoin(bt_index_dir,'{seq_index}.{num}.ebwt'), seq_index = seq_index_list, num = [1,2,3,4]),
		in_rev_indexes = expand(pjoin(bt_index_dir,'{seq_index}.rev.{num}.ebwt'), seq_index = seq_index_list, num = [1,2])
	output:
		out_seq_map = expand(pjoin(bowtie_out_dir, '{lib}_{seq_index}_mapped.txt'), lib=lib_list, seq_index = seq_index_list)
	run:
		for lib in lib_list:
			mapping_list = seq_index_list.copy() #reset mapping_list variable as seq_index_list at each lib
			print(mapping_list)
			current_lib = pjoin(lib_dir, lib + '_final_col_ntfilter.fasta') #initialize library variable
			while mapping_list: #while index_list is not empty enter sequential mapping 
				current_index = mapping_list.pop(0) #while mapping_list is not empty
				unmapped_fasta = pjoin(bowtie_out_dir, lib + '_' + current_index + '_unmapped.fasta')
				shell('bowtie -f --un {unmapped_fasta} -v 0 -p 1 -t {bt_index_dir}/{current_index} \
				{current_lib} > {bowtie_out_dir}/{lib}_{current_index}_mapped.txt')
				if (current_lib != pjoin(lib_dir, lib + '_final_col_ntfilter.fasta')):
					shell('rm {current_lib}')
				current_lib = unmapped_fasta #reset current_lib variable


rule Generate_rawcounts:
	input:
		in_seq_map = expand(pjoin(bowtie_out_dir, '{lib}_{index}_mapped.txt'),lib=lib_list, index = full_index_list),
	output:
		out_rawcounts = pjoin(counts_out_dir, '{index}_rawcounts.txt')
	run:
		if ((lengthfilter == 'true')  & (wildcards.index in siRNA_index_list)):
			GatherLengthfilterBowtieCounts(lib_list,
				wildcards.index,
				bowtie_out_dir,
				counts_out_dir,
				minlength,
				maxlength)
		else:
			GatherBowtieCounts(lib_list,
				wildcards.index,
				bowtie_out_dir,
				counts_out_dir)




'''rule Generate_main_figures:
	input:
		expand(pjoin(counts_out_dir, '{index}_rawcounts.txt'), index = full_index_list)
	output:
		out_mainfigs = pjoin(figures_out_dir,'{index_type}_{fig_type}_dataset={dataset}_runID={runID}_.png'),
	run:
		shell('Rscript --vanilla Scripts/R/main_figures.R {dataset} {runID}') 


rule Generate_maplog_figure:
	input:
		expand(pjoin(counts_out_dir, '{index}_rawcounts.txt'), index = full_index_list)
	output:
		pjoin(maplog_out_dir, dataset + '_maplog.png')
	run:
		shell('Rscript --vanilla Scripts/R/main_figures.R {dataset} {runID}') 

'''

