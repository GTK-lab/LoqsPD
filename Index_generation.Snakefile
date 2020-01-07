'''
Pre-conditions:
input_dir = "Input/Index_Generation/"
1)Folders that need to be created/ stored with necessary files
	a) /Downloads -> With: Gff file, Genome fasta , 4 citation bedfiles, encode cis-nat bedfile, RM output
	Chr1dmerdna.fasta, spikeinset01.fasta, Genome

	b) /Configfiles -> With: main_config.yaml

2) Usage
	snakemake -s index_generation.Snakefile 

'''

from os.path import join as pjoin
from Scripts.Python.helpers import ExtractExonWithGeneID,ExtractFtFromDB,ExtractGeneIDFromGff, Gff2BedMod, CisNatBedMod, RmBedMod
import pandas as pd
import gffutils 
import sys

configfile: srcdir("Configfiles/main_config.yaml")
downloads_dir = 'Input/Downloads'
extracted_gffbed_dir = 'Output/Index_generation/Extracted_gff_bed'
final_bed_dir = 'Output/Index_generation/Final_bedfiles'
final_fasta_dir = 'Output/Index_generation/Final_fastafiles'
bt_index_dir = 'Output/Index_generation/Bowtie_indexes'

genome_vers = config['genome_vers']
index_list = list(config['index_dict'].keys())
#rm_ft_list= config['rm_ft_list']
citation_ft_list= config['citation_ft_list']
downloaded_fasta_list = config['downloaded_fasta_list']
genome_path = pjoin(downloads_dir, genome_vers + '.fasta')


rule all:
	input:
		#expand(pjoin(final_fasta_dir,'{rm_ft}RM.fasta'), rm_ft = rm_ft_list),
		expand(pjoin(final_fasta_dir, '{citation_ft}.fasta'), citation_ft =citation_ft_list),
		expand(pjoin(final_fasta_dir, '{downloaded}.fasta'), downloaded = downloaded_fasta_list),  
		expand(pjoin(bt_index_dir,'{index}.{num}.ebwt'), index= index_list, num = [1,2,3,4]),
		expand(pjoin(bt_index_dir,'{index}.rev.{num}.ebwt'), index= index_list, num = [1,2])

#assume genomeDB created 
rule Extractft_geneID_gff:
	input:
		pjoin(extracted_gffbed_dir, genome_vers + '.db')
	output:
		out_gff = pjoin(extracted_gffbed_dir,'{ft}.gff'),
		out_geneID = pjoin(extracted_gffbed_dir,'{ft}_geneID.txt'),
		out_exon_gff = pjoin(extracted_gffbed_dir,'{ft}_exon.gff'),
	run:
		db = gffutils.FeatureDB(str(input), keep_order=True)
		ExtractFtFromDB(db, wildcards.ft, output.out_gff)
		ExtractGeneIDFromGff(output.out_gff, output.out_geneID)
		ExtractExonWithGeneID(db, 
							  output.out_geneID,
							  wildcards.ft, 
							  output.out_exon_gff)

rule Extract_hpRNA_features:
	input:
		pjoin(extracted_gffbed_dir, genome_vers + '.db')
	output:
		out_gene_gff = pjoin(extracted_gffbed_dir, 'gene.gff'),
		out_hpRNA_gff = pjoin(extracted_gffbed_dir, 'hpRNA.gff')
	run:
		db = gffutils.FeatureDB(str(input), keep_order=True)
		ExtractFtFromDB(db, 'gene', output.out_gene_gff)
		shell('grep "Name=hpRNA" {output.out_gene_gff} > {output.out_hpRNA_gff}')


rule Gff2bed_and_merge_for_exon_TE:
	input: 
		in_gff = pjoin(extracted_gffbed_dir,'{bed_ft}.gff')
	output:
		out_bed = pjoin(extracted_gffbed_dir,'{bed_ft}.bed'),
		out_merge_bed = pjoin(final_bed_dir, '{bed_ft}_merge.bed'),
		out_fasta = pjoin(final_fasta_dir,'{bed_ft}_merge.fasta')
	run:
		shell("gff2bed < {input.in_gff} > {output.out_bed}")
		Gff2BedMod(output.out_bed)
		shell('bedtools merge -s -c 4,5,6 -o distinct -i {output.out_bed} > {output.out_merge_bed}')
		shell('bedtools getfasta -name -fi {genome_path} -bed {output.out_merge_bed} -fo {output.out_fasta}')


rule Extract_mRNA_exon_stranded_features: 
	input: 
		pjoin(final_bed_dir,'mRNA_exon_merge.bed')
	output:
		out_stranded_mRNA_exon = pjoin(final_bed_dir, 'mRNA_exon_merge_{strand}.bed')
	run: 
		exon_bed = pd.read_csv(str(input), sep='\t', header=None, index_col=None)
		strand_sign = '+' if wildcards.strand == 'plus' else '-'
		strand_features = exon_bed.loc[exon_bed[5] == strand_sign]
		strand_features.to_csv(output.out_stranded_mRNA_exon, sep='\t', header=None, index=False)


rule Generate_antiexon_mRNA_mRNA_bed_fasta:
	input:
		mRNA_exon_plus = pjoin(final_bed_dir, 'mRNA_exon_merge_plus.bed'),
		mRNA_exon_minus = pjoin(final_bed_dir, 'mRNA_exon_merge_minus.bed'),
	output:
		 out_bed = pjoin(final_bed_dir, 'antiexon_mRNA_mRNA.bed'),
		 out_fasta = pjoin(final_fasta_dir, 'antiexon_mRNA_mRNA.fasta')
	run:
		shell('bedtools intersect -S -wb -a {input.mRNA_exon_plus} -b {input.mRNA_exon_minus} > {output.out_bed}')
		CisNatBedMod(str(output.out_bed))
		shell('bedtools getfasta -name -fi {genome_path} -bed {output.out_bed} -fo {output.out_fasta}')


rule Generate_mRNA_other_cisnat_bed_fasta:
	input:
		in_other_bed = pjoin(final_bed_dir, '{antiexon_ft}_exon_merge.bed'),
		in_mRNA_bed  = pjoin(final_bed_dir, 'mRNA_exon_merge.bed'),
	output:
		out_bed = pjoin(final_bed_dir, 'antiexon_mRNA_{antiexon_ft}.bed'),
		out_fasta = pjoin(final_fasta_dir, 'antiexon_mRNA_{antiexon_ft}.fasta')
	run:
		shell('bedtools intersect -S -wb -a {input.in_mRNA_bed} -b {input.in_other_bed} > {output.out_bed}')
		CisNatBedMod(str(output.out_bed))
		shell('bedtools getfasta -name -fi {genome_path} -bed {output.out_bed} -fo {output.out_fasta}')


rule Generate_RM_bed_fasta: 
	input:
		pjoin(downloads_dir, genome_vers + '.fasta.RMout')
	output:
		out_bed = pjoin(final_bed_dir,'{rm_ft}RM.bed'),
		out_fasta = pjoin(final_fasta_dir,'{rm_ft}RM.fasta')
	run:
		RmBedMod(str(input), output.out_bed, wildcards.rm_ft)
		shell('bedtools getfasta -name -fi {genome_path} -bed {output.out_bed} -fo {output.out_fasta}')


rule Generate_citation_fasta:
	input:
		in_bed = expand(pjoin(downloads_dir, '{citation_ft}.bed'), citation_ft =citation_ft_list)
	output:
		out_fasta = expand(pjoin(final_fasta_dir, '{citation_ft}.fasta'), citation_ft =citation_ft_list)
	run:
		for citation_lib in citation_ft_list:
			shell('cp {downloads_dir}/{citation_lib}.bed {final_bed_dir}') #sneaky line added to copy over bedfiles
			shell('bedtools getfasta -name -fi {genome_path} -bed {final_bed_dir}/{citation_lib}.bed -fo {final_fasta_dir}/{citation_lib}.fasta')


rule Move_downloaded_fasta_files:
	input:
		in_fasta =expand(pjoin(downloads_dir, '{downloaded}.fasta'), downloaded = downloaded_fasta_list)
	output:
		out_fasta = expand(pjoin(final_fasta_dir, '{downloaded}.fasta'), downloaded = downloaded_fasta_list)
	run: 
		shell('cp {input.in_fasta} {final_fasta_dir}')


rule Generate_bt_indexes:
	input:
		expand(pjoin(final_fasta_dir, '{index}.fasta'), index= index_list)
	output:
		expand(pjoin(bt_index_dir,'{index}.{num}.ebwt'), index= index_list, num = [1,2,3,4]),
		expand(pjoin(bt_index_dir,'{index}.rev.{num}.ebwt'), index= index_list, num = [1,2])
	run:
		for index in index_list:
			in_fasta_path = pjoin(final_fasta_dir, index + '.fasta')
			out_index_path = pjoin(bt_index_dir, index )
			shell('bowtie-build {in_fasta_path} {out_index_path}')


