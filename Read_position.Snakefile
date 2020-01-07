'''
Run as:

snakemake -s Read_position.Snakefile --config dataset=''  runID ='' lengthfilter ='' minlength = '' maxlength= ''

- 'dataset' and 'runID': Strings which corresponds to the directory where the bowtie output is stored, i.e `Output/dataset/Mapping/runID=runID/Bowtie_output` 
- 'lengthfilter': Boolean true/false that determines whether length filtering is performed 
- 'minlength': Integer which sets the minimum (greater than or equal to) length for a mapped read to be recognized as a count. 
- 'minlength': Integer which sets the maximum (less than or equal to) length for a mapped read to be recognized as a count. 


**Note that target bowtie output are set within the script, under FC_list
'''

import os
from pathlib import Path
from os.path import join as pjoin
from Scripts.Python.helpers import MakeGenomeMapLog, SortMapFiles
import pandas as pd
import sys

configfile: srcdir("Configfiles/main_config.yaml")

dataset = config['dataset']
runID=  config['runID']
lengthfilter = config['lengthfilter']

if lengthfilter == 'true' :
	minlength = config['minlength']
	maxlength = config['maxlength']
	bowtie_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/Bowtie_output'
	filter_out_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/' + str(minlength) + '-' + str(maxlength) + 'filter_Bowtie_output'
	strand_out_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/' + str(minlength) + '-' + str(maxlength) + 'strand_output'
else:
	bowtie_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/Bowtie_output'
	filter_out_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/Bowtie_output'
	strand_out_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/' + 'strand_output' 

other_counts_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/Counts_output'



siRNA_list = ['hpRNA_merge','transposable_element','encodecisnat',
'01_okamura', '02_czech', '03_kawamura', '04_ghildiyal',
'antiexon_mRNA_mRNA','antiexon_mRNA_pseudogene','antiexon_mRNA_ncRNA','antiexon_mRNA_snoRNA', 'antiexon_mRNA_tRNA']
sRNA_list = siRNA_list + ['miRNA_merge']
ft_list = ['hpRNA_merge', 'encodecisnat']
misc_list = ['dmel-all-chromosome-r6.22','spikeinset01']
lib_list = list(config[dataset + '_lib_dict'].keys())



rule all:
	input:
		expand(pjoin(filter_out_dir,'{lib}_{feature}_mapped.txt'), feature= ft_list, lib = lib_list),
		expand(pjoin(strand_out_dir,'{feature}','{lib}_{feature}_{strand}_readpos.txt'), feature= ft_list, lib = lib_list, strand = ["plus","minus"])

rule lengthfilter: 
	input:
		pjoin(bowtie_dir,'{lib}_{feature}_mapped.txt')
	output:
		pjoin(filter_out_dir,'{lib}_{feature}_mapped.txt')
	run:
		in_file = str(input)
		if os.stat(in_file).st_size != 0:
			mapped = pd.read_csv(in_file,sep='\t', header=None, index_col=None)
			mapped.length = mapped[0].str.split('_', expand = True)[3].astype('str')
			mapped.length = pd.to_numeric(mapped.length.str.split('nt', expand=True)[0])
			mapped_filt = mapped[(mapped.length >= minlength) & (mapped.length <= maxlength)]
			mapped_filt.to_csv(str(output), sep='\t', index=False, header=None)
		else:
			empty_df = pd.DataFrame()
			empty_df.to_csv(str(output), sep='\t', index=False, header= None)


rule split_strand: 
	input:
		pjoin(filter_out_dir,'{lib}_{feature}_mapped.txt')
	output:
		plus_pos = pjoin(strand_out_dir,'{feature}','{lib}_{feature}_plus_readpos.txt'),
		minus_pos = pjoin(strand_out_dir,'{feature}','{lib}_{feature}_minus_readpos.txt')

	run:
		file = str(input)
		if os.stat(file).st_size != 0:
			mapped = pd.read_csv(file,sep='\t', header=None, index_col=None)
			mapped['counts'] = mapped[0].str.split('_', expand = True)[2].astype('float')
			mapped['length'] = pd.to_numeric(mapped[0].str.split('_', expand = True)[3].str.split('nt', expand = True)[0])
			minus = mapped[mapped[1].str.contains("-")][["counts","length",3,2]]
			plus = mapped[~mapped[1].str.contains("-")][["counts","length",3,2]]
			minus.columns = ["counts","length","sense_pos","loci"]
			plus.columns = ["counts","length","sense_pos","loci"]
			plus_counts = plus.groupby(['loci',"length","sense_pos"])["counts"].sum().reset_index(name='group_counts')
			minus_counts = minus.groupby(['loci',"length","sense_pos"])["counts"].sum().reset_index(name='group_counts')
			plus_counts['lib'] = wildcards.lib
			minus_counts['lib'] = wildcards.lib

			plus_counts.to_csv(str(output.plus_pos), sep='\t', index=False)
			minus_counts.to_csv(str(output.minus_pos), sep='\t', index=False)
		else: 
			empty_df = pd.DataFrame()
			empty_df.to_csv(str(output.plus_pos), sep='\t', index=False)
			empty_df.to_csv(str(output.minus_pos), sep='\t', index=False)






