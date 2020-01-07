'''
Run as:
snakemake -s Externallib_preprocessing.Snakefile --config dataset=''

- 'dataset' corresponds to the name of directory where the raw libraries are stored i.e `Input/Libraries/dataset`

'''
import sys
import pandas as pd 
import os 
from os.path import join as pjoin

configfile: srcdir("Configfiles/main_config.yaml")
dataset = config['dataset']

lib_list = list(config[dataset + '_lib_dict'].keys())

raw_libdir = 'Input/Libraries/'+ dataset 
final_libdir = 'Output/'+ dataset + '/Preprocessed_libs/Final_preprocessed_libs'


rule all:
	input: 
		expand(pjoin(final_libdir, '{lib}_final_col_filter.fasta'), lib = lib_list),
		expand(pjoin(final_libdir,'{lib}_final_col_ntfilter.fasta'), lib=lib_list)

rule GSE17171:
	input:
		expand(pjoin(raw_libdir, '{lib_condition}_{maptype}.csv'), 
			lib_condition = ['GSM430035_lacZ','GSM430034_r2d2','GSM430033_loqs-ORF','GSM430032_dcr-2','GSM430031_dcr-1','GSM430030_mock'],
			maptype = ['mappers', 'non_mappers'])
	output:
		pjoin('Output/GSE17171/Preprocessed_libs/Final_preprocessed_libs', '{lib}_final_col_filter.fasta'),

	run:
		file_lst = [file for file in os.listdir(raw_libdir) if wildcards.lib in file]
		mapped = pd.read_csv(pjoin(raw_libdir, file_lst[0]))
		unmapped = pd.read_csv(pjoin(raw_libdir, file_lst[1]))
		all_reads =  pd.concat([mapped,unmapped])
		sorted_reads = all_reads.sort_values(by=['READS'], ascending=False)
		sorted_reads['LENGTH'] = sorted_reads['SEQ'].str.len()
		filtered_reads = sorted_reads.loc[(sorted_reads['LENGTH'] >= 18) & (sorted_reads['LENGTH'] <= 30)]
		filtered_reads = filtered_reads.reset_index(drop=True)
		filtered_reads['RANK'] = filtered_reads.index + 1
		filtered_reads['READS'] = '>' + filtered_reads['RANK'].astype(str) + '_' + wildcards.lib + '_' + filtered_reads['READS'].astype(str)

		df_list = []
		for row in range(len(filtered_reads)):
				df_list.append(filtered_reads.loc[row]['READS'])
				df_list.append(filtered_reads.loc[row]['SEQ'])
		fasta_df = pd.DataFrame(df_list)
		fasta_df.to_csv('temp.fasta', sep='\t', index=False, header = False)
		shell('fastx_artifacts_filter -v -i temp.fasta -o {output}') 
		shell('rm temp.fasta')



#Note, renamed  GSM643930_r2d2_21.txt to GSM643930_r2d2_19_25.txt i cant take this is shit
rule GSE26230:
	input:
		expand(pjoin(raw_libdir, '{lib_condition}_19_25.txt'), 
			lib_condition = ['GSM643927_loqsCyo','GSM643928_loqs','GSM643929_r2d2Cyo','GSM643930_r2d2'])
	output:
		pjoin('Output/GSE26230/Preprocessed_libs/Final_preprocessed_libs', '{lib}_final_col_filter.fasta')

	run:
		for file in os.listdir(raw_libdir):
			if file.startswith(wildcards.lib) and file.endswith('_19_25.txt'):
				fastq = pd.read_csv(pjoin(raw_libdir,file),header=None, index_col=None) 
				fasta_tab = pd.DataFrame(columns=['READS','ID'])
				fasta_tab['READS'] = fasta_tab['READS'] = fastq.loc[1::4].reset_index(drop=True)[0]
				fasta_tab['ID'] = '>' + wildcards.lib
				fasta_tab = fasta_tab.loc[~fasta_tab['READS'].str.contains('N')]
				fasta_tab = fasta_tab.loc[~fasta_tab['READS'].str.contains('a|t|c|g')].reset_index(drop=True)
				
				df_list = []
				for row in range(len(fasta_tab)):
					df_list.append(fasta_tab.loc[row]['ID'])
					df_list.append(fasta_tab.loc[row]['READS'])

				fasta_df = pd.DataFrame(df_list)
				fasta_df.to_csv('raw.fasta', sep='\t', index=False, header = False)
				shell('fastx_collapser -v -i raw.fasta | fastx_artifacts_filter -v -o {output}')

				col_fasta = pd.read_csv(str(output),  header = None, index_col=None)

				for row in range(len(col_fasta)):
					data = col_fasta.loc[row]
					if data.str.startswith('>')[0]:
						counts = data.str.split('-')[0][1]
						col_fasta.loc[row] = '>'+ str(row+1) + '_' + wildcards.lib + '_' +  str(counts)
				col_fasta.to_csv(str(output),  header = None, index=False)

		
rule GSE37443Head:
	input:
		expand(pjoin(raw_libdir, '{lib_condition}.txt'), 
			lib_condition = ['GSM919400_wIR++Un', 'GSM919401_wIRloqsKOPAUn', 'GSM919402_wIRloqsKOPBUn', 'GSM919403_wIRloqsKOPABUn', 
				'GSM919404_wIRloqsKOPADUn', 'GSM919405_wIRloqsKOPBDUn', 'GSM919406_wIRloqsKOPABDUn', 'GSM919407_wIRloqsKOCyOUn',
				'GSM919408_wIRloqsHypUn', 'GSM919409_wIRloqsHypCyOUn'])
	output:
		pjoin('Output/GSE37443Head/Preprocessed_libs/Final_preprocessed_libs', '{lib}_final_col_filter.fasta')

	run:
		file_lst = [file for file in os.listdir(raw_libdir) if wildcards.lib in file]
		fasta_tab = pd.read_csv(pjoin(raw_libdir, file_lst[0]), sep='\t', names=['READS', 'COUNTS'], index_col=None)
		fasta_tab['READS'] = fasta_tab['READS'].str.replace('.','N')
		fasta_tab = fasta_tab.loc[~fasta_tab['READS'].str.contains('N')]
		fasta_tab = fasta_tab.loc[~fasta_tab['READS'].str.contains('a|t|c|g')]
		fasta_tab = fasta_tab.sort_values(by='COUNTS', ascending=0).reset_index(drop=True)
		fasta_tab['RANK'] = fasta_tab.index + 1
		fasta_tab['COUNTS'] = '>' + fasta_tab['RANK'].astype(str) + '_' + wildcards.lib + '_' + fasta_tab['COUNTS'].astype(str)

		df_list = []
		for row in range(len(fasta_tab)):
			df_list.append(fasta_tab.loc[row]['COUNTS'])
			df_list.append(fasta_tab.loc[row]['READS'])
		fasta_df = pd.DataFrame(df_list)
		temp_fasta = pjoin(raw_libdir,'temp.fasta')
		fasta_df.to_csv(temp_fasta,sep='\t', header=None, index=False)

		shell('fastx_clipper -v -a "TCGTATGCCGTCTTCTGCTTG" -l 18 -c -i {temp_fasta} | fastx_artifacts_filter -v -o {output}')

rule GSE37443Ovary:
	input:
		expand(pjoin(raw_libdir, '{lib_condition}.txt'), 
			lib_condition = ['GSM919410_OVw1118Un', 'GSM919411_OVloqsKOPAUn',
				'GSM919412_OVloqsKOPBUn', 'GSM919413_OVloqsKOPABUn', 'GSM919414_OVloqsKOPADUn', 'GSM919415_OVloqsKOPBDUn',
				'GSM919416_OVloqsKOPABDUn', 'GSM919417_OVloqsKOCyOUn', 'GSM919418_OVloqsHypUn', 'GSM919419_OVloqsHypCyOUn'])
	output:
		pjoin('Output/GSE37443Ovary/Preprocessed_libs/Final_preprocessed_libs', '{lib}_final_col_filter.fasta')

	run:
		for file in os.listdir(raw_libdir):
			if file.startswith(wildcards.lib):
				print(file)
				fasta_tab = pd.read_csv(pjoin(raw_libdir, file), sep='\t', names=['READS', 'COUNTS'], index_col=None)
				fasta_tab['READS'] = fasta_tab['READS'].str.replace('.','N')
				fasta_tab = fasta_tab.loc[~fasta_tab['READS'].str.contains('N')]
				fasta_tab = fasta_tab.loc[~fasta_tab['READS'].str.contains('a|t|c|g')]
				fasta_tab = fasta_tab.sort_values(by='COUNTS', ascending=0).reset_index(drop=True)
				fasta_tab['RANK'] = fasta_tab.index + 1
				fasta_tab['COUNTS'] = '>' + fasta_tab['RANK'].astype(str) + '_' + wildcards.lib + '_' + fasta_tab['COUNTS'].astype(str)
				print('check1')

				df_list = []
				for row in range(len(fasta_tab)):
					df_list.append(fasta_tab.loc[row]['COUNTS'])
					df_list.append(fasta_tab.loc[row]['READS'])
				fasta_df = pd.DataFrame(df_list)
				temp_fasta = pjoin(raw_libdir, wildcards.lib + '_temp.fasta')
				print('check3')
				fasta_df.to_csv(temp_fasta,sep='\t', header=None, index=False)
				shell('fastx_clipper -v -a "TCGTATGCCGTCTTCTGCTTG" -l 18 -c -i {temp_fasta} | fastx_artifacts_filter -v -o {output}')


rule Adding_read_length:
	input:
		pjoin(final_libdir, '{lib}_final_col_filter.fasta')

	output:
		pjoin(final_libdir, '{lib}_final_col_ntfilter.fasta')

	run:
		temp_tab = pjoin(final_libdir, wildcards.lib +'_temp.tab')

		shell('cat {input} | paste - - > {temp_tab}')
		fasta_tab = pd.read_csv(temp_tab, sep='\t', names=['ID', 'READ'], index_col=None)
		fasta_tab['LENGTH'] = fasta_tab['READ'].str.len()
		fasta_tab['ID'] = fasta_tab['ID'].astype(str) + '_' + fasta_tab['LENGTH'].astype(str) + 'nt'

		df_list = []
		for row in range(len(fasta_tab)):
			df_list.append(fasta_tab.loc[row]['ID'])
			df_list.append(fasta_tab.loc[row]['READ'])
		fasta_df = pd.DataFrame(df_list)
		fasta_df.to_csv(str(output),sep='\t', header=None, index=False)
		shell('rm {temp_tab}')
