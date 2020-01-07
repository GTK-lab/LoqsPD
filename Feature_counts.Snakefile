'''
Run as:

snakemake -s Feature_counts.Snakefile --config dataset='' runID='' lengthfilter='' minlength='' maxlength=''

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

#Defining --config variables
dataset = config['dataset']
runID  = config['runID'] #default named as Liling
lengthfilter = config['lengthfilter']


if lengthfilter == 'true' :
	minlength = config['minlength']
	maxlength = config['maxlength']


lib_list = list(config[dataset + '_lib_dict'].keys())
FC_list = ['hpRNA_merge','transposable_element']
## ['transposable_element']
## ['antiexon_mRNA_mRNA', 'antiexon_mRNA_tRNA', 'antiexon_mRNA_pseudogene', 'antiexon_mRNA_snoRNA','antiexon_mRNA_ncRNA']
##['antiexon_mRNA_mRNA', 'antiexon_mRNA_tRNA', 'antiexon_mRNA_pseudogene', 'antiexon_mRNA_snoRNA','antiexon_mRNA_ncRNA']


bowtie_output_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/Bowtie_output'
FC_out_dir = 'Output/' + dataset + '/Mapping/runID=' + runID + '/FC_output'


rule all:
	input:
		#expand(pjoin('Flags', dataset + '_runID=' + runID + '_{feature}_FC_dirs_initialized.done'), feature = FC_list),
		expand(pjoin(FC_out_dir,'{feature}_FC','{lib}_{feature}_locicounts.txt'), feature= FC_list, lib = lib_list)

'''
#due to me seemingly not being able to generate more than 2 conditoinal outputs in Intialize directories
#We are gonna have to use this i guess. 
rule Make_FC_dirs:
	output:
		touch(pjoin('Flags', dataset + '_runID=' + runID + '_{feature}_FC_dirs_initialized.done'))
	run:
		dirs = ['Output/' + dataset + '/Mapping/runID=' + runID + '/FC_output/' + wildcards.feature + '_FC',
				'Output/' + dataset + '/Mapping/runID=' + runID + '/Figures_output/' + wildcards.feature + '_FC'
				]
		for directory in dirs:
			if not Path(directory).is_dir():
				os.makedirs(directory)
'''

rule counts_by_loci:
	input:
		pjoin(bowtie_output_dir,'{lib}_{feature}_mapped.txt')
	output:
			pjoin(FC_out_dir,'{feature}_FC','{lib}_{feature}_locicounts.txt') 
	run:
		file = str(input)
		if os.stat(file).st_size != 0:
			mapped = pd.read_csv(file,sep='\t', header=None, index_col=None)
			if lengthfilter == 'true':
				mapped.length = mapped[0].str.split('_', expand = True)[3].astype('str')
				mapped.length = pd.to_numeric(mapped.length.str.split('nt', expand=True)[0])
				mapped= mapped[(mapped.length >= minlength) & (mapped.length <= maxlength)]
			counts = mapped[0].str.split('_', expand = True)[2].astype('int')
			mapped[0] = counts 
			df = mapped.groupby([2], as_index=False)[0].sum()
			df[3] = 'Counts'
			pivotdf = df.pivot(columns=2, index=3, values = 0)
			pivotdf.to_csv(str(output), sep='\t', index=False, header= list(pivotdf))
		else:
			empty_df = pd.DataFrame()
			empty_df.to_csv(str(output), sep='\t', index=False, header= None)