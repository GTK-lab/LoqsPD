'''
Run as:

snakemake -s Multi_sequential_mapping.Snakefile --config dataset=''  runID ='' minlength = '' maxlength= ''


- 'dataset': String which corresponds to the directory where the pre-processed libraries are stored, i.e `Output/dataset/Preprocessed_libs/Final_preprocessed_libs`
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
dataset = config['dataset']
runID  = config['runID']
lengthfilter = config['lengthfilter']


if lengthfilter == 'true' :
	minlength = config['minlength']
	maxlength = config['maxlength']


#Setting index variables and directories
bt_index_dir = 'Output/Index_generation/Bowtie_indexes'
#seq_index_list = ['dmel-all-chromosome-r6.22'] + config[runID + '_index_list']
seq_index_list = config[runID + '_index_list']
nonseq_index_list = config['nonseq_index_list']
full_index_list = seq_index_list + nonseq_index_list
foo_index_list = seq_index_list + ["spikeinset01"]
siRNA_index_list = config['siRNA_index_list'] 
index_dict = config['index_dict']
index_type_list = list(set(index_dict.get(x) for x in seq_index_list))
index_type_list.append('Combined')

#Genome version 
genome_vers = config['genome_vers']

#Setting lib_list based on {dataset}
lib_list = list(config[dataset + '_lib_dict'].keys())

lib_dir = 'Output/' + dataset + '/Preprocessed_libs/Final_preprocessed_libs'
bowtie_out_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/Bowtie_output'
counts_out_dir = 'Output/' + dataset + '/Mapping/runID=' +  runID + '/Counts_output'


rule all:
	input:
		#expand(pjoin(maplog_out_dir, '{lib}_maplog.txt'), lib = lib_list),
		expand(pjoin(counts_out_dir, '{index}_rawcounts.txt'), index = foo_index_list)
		#expand(pjoin(bowtie_out_dir,'{lib}_' + genome_vers + '_unmapped.fasta'), lib=lib_list),
		#pjoin(maplog_out_dir, dataset + '_maplog.png'),


rule multi_equential_mapping:
	input:
		in_libs = expand(pjoin(lib_dir,'{lib}_final_col_ntfilter.fasta'), lib=lib_list),
		in_fw_indexes = expand(pjoin(bt_index_dir,'{seq_index}.{num}.ebwt'), seq_index = seq_index_list, num = [1,2,3,4]),
		in_rev_indexes = expand(pjoin(bt_index_dir,'{seq_index}.rev.{num}.ebwt'), seq_index = seq_index_list, num = [1,2])
	output:
		out_seq_map = expand(pjoin(bowtie_out_dir, '{lib}_{seq_index}_mapped.txt'), lib=lib_list, seq_index = seq_index_list)
	run:
		for lib in lib_list:
			print(bowtie_out_dir)
			mapping_list = seq_index_list.copy() #reset mapping_list variable as seq_index_list at each lib
			print(mapping_list)
			current_lib = pjoin(lib_dir, lib + '_final_col_ntfilter.fasta') #initialize library variable
			while mapping_list: #while index_list is not empty enter sequential mapping 
				current_index = mapping_list.pop(0) #while mapping_list is not empty
				unmapped_fasta = pjoin(bowtie_out_dir, lib + '_' + current_index + '_unmapped.fasta')
				shell('bowtie -f --un {unmapped_fasta} --all -v 0 -p 1 -t {bt_index_dir}/{current_index} \
				{current_lib} > {bowtie_out_dir}/{lib}_{current_index}_mapped.txt')
				current_file = pjoin(bowtie_out_dir, lib + '_'+ current_index + '_mapped.txt')
				if os.stat(current_file).st_size != 0:
					mapped = pd.read_csv(current_file,sep='\t', header=None, index_col=None)
					counts = mapped[0].str.split('_', expand = True)[2].astype('int') / (mapped[6].astype('int') + 1)
					header =  mapped[0].str.split('_', expand = True)
					mapped[0] = header[0].astype(str) + '_' + header[1].astype(str) + '_' + counts.astype(str)+ '_' + header[3].astype(str)
					mapped.to_csv(current_file, sep='\t', header=None, index=False)
				if (current_lib != pjoin(lib_dir, lib + '_final_col_ntfilter.fasta')):
					shell('rm {current_lib}')
				if current_index != 'dmel-all-chromosome-r6.22':
					current_lib = unmapped_fasta #reset current_lib variable
 
rule Spikein_mapping:
	input:
		in_libs = pjoin(lib_dir,'{lib}_final_col_ntfilter.fasta'),
		in_fw_indexes = expand(pjoin(bt_index_dir,'spikeinset01.{num}.ebwt'), num = [1,2,3,4]),
		in_rev_indexes = expand(pjoin(bt_index_dir,'spikeinset01.rev.{num}.ebwt'), num = [1,2])
	output:
		multi_mapped = pjoin(bowtie_out_dir,'{lib}_spikeinset01_mapped.txt'),
	run:
		index_path = pjoin(bt_index_dir, 'spikeinset01')
		shell('bowtie -f --all -v 0 -p 1 -t {index_path} {input.in_libs} > {output.multi_mapped}')


rule Generate_rawcounts:
	input:
		in_seq_map = expand(pjoin(bowtie_out_dir, '{lib}_{index}_mapped.txt'),lib=lib_list, index = foo_index_list),
	output:
		out_rawcounts = pjoin(counts_out_dir, '{index}_rawcounts.txt')
	run:
		if ( (lengthfilter == 'true') & (wildcards.index in siRNA_index_list)):
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
		