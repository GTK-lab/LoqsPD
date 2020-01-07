'''
Run as:

snakemake -s Okamuralib_preprocessing.Snakefile

'''
from os.path import join as pjoin
import pandas as pd

configfile: srcdir('Configfiles/main_config.yaml')
raw_libdir = 'Input/Libraries/Okamura'
inter_libdir = 'Output/Okamura/Preprocessed_libs/Intermediate_libs'
final_libdir = 'Output/Okamura/Preprocessed_libs/Final_preprocessed_libs'

lib_list = list(config['Okamura_lib_dict'].keys())

adapter = config['adapter']

rule all:
	input:
		expand(pjoin(final_libdir,'{lib}_final_col_ntfilter.fasta'), lib=lib_list)

##Takes around 10 mins to run
#sed line is very handy 
rule Clip_convert_to_fasta_and_label_lib :
	input:
		in_fastq = pjoin(raw_libdir,'{lib}.fastq')
	output:
		out_fasta = pjoin(inter_libdir,'{lib}_clip.fasta')
	run:
		shell('fastx_clipper -a {adapter} -c -v -i {input.in_fastq} | fastq_to_fasta -v -r | \
				sed "s/^>/>{wildcards.lib}_/" > {output.out_fasta}')

#After collapsing, reads sorted via counts.
#>1_B133_2416647
#TGAGATCATTTTGAAAGCTGATT
rule Collapse_repeats:
	input:
		in_fasta = pjoin(inter_libdir,'{lib}_clip.fasta')
	output:
		out_fasta = pjoin(inter_libdir,'{lib}_clip_collapse.fasta')
	run:
		shell("fastx_collapser -v -i {input.in_fasta} -o {output.out_fasta}")


rule Convert_to_tab_and_filter_length:
	input:
		in_fasta = pjoin(inter_libdir,'{lib}_clip_collapse.fasta')
	output:
		out_tab =pjoin(inter_libdir,'{lib}_len.tab')
	run:
		shell('fasta_formatter -t -i {input.in_fasta}> {output.out_tab}')
		tab_data = pd.read_csv(output.out_tab, sep='\t',header=None, index_col=None)
		tab_data['length'] = tab_data[1].str.len()
		filtered_data = tab_data.query('18 <= length <= 30')
		filtered_data.to_csv(output.out_tab, sep='\t', header=None, index=False)

#using iterrows and indexing the row item for identifier + read 
#Note as of 02/01/2019, the last field is now the length the read
#read_ID = rank_lib_counts_length 
rule Convert_to_fasta:
	input:
		in_tab = pjoin(inter_libdir,'{lib}_len.tab')
	output:
		out_fasta = pjoin(inter_libdir,'{lib}_len.fasta')
	run:
		sys.stdout = open(output.out_fasta, 'w')
		tab_data = pd.read_csv(input.in_tab, sep = '\t', header = None, index_col= None)
		count_fields = tab_data[0].str.split('-',expand=True)
		tab_data[0] = '>'+ count_fields[0].astype(str) + '_' + wildcards.lib + '_' +  count_fields[1].astype(str) + '_' + tab_data[2].astype(str) + 'nt'
		for index,row in tab_data.iterrows():
			data = row[:2]
			for line in data:
				print(line)
		sys.stdout.close()


rule Artifact_filter:
	input:
		in_fasta = pjoin(inter_libdir,'{lib}_len.fasta')
	output:
		out_fasta = pjoin(final_libdir,'{lib}_final_col_ntfilter.fasta')
	run:
		shell('fastx_artifacts_filter -v -i {input.in_fasta} > {output.out_fasta}')
		shell('echo "Final read identiier count:" grep -c {wildcards.lib} {output.out_fasta}')


