import os 
from os.path import join as pjoin
import pandas as pd
import gffutils
import sys



def ExtractFtFromDB(db, ft,out_gff):
	sys.stdout = open(out_gff,'w')
	for feature in db.features_of_type(ft):
		print(feature)
	sys.stdout.close()

def ExtractGeneIDFromGff(gff,out_geneID):
	gffdata= pd.read_csv(gff,sep='\t', header=None, index_col=None)
	geneID_field = (gffdata[8].str.split(';',expand=True))[3]
	geneID = (geneID_field.str.split('=',expand=True))[1].drop_duplicates()
	geneID.to_csv(out_geneID,sep='\t', header=None, index=False)

def ExtractExonWithGeneID(db,geneID,ft,out_gff):
	geneIDs = [line.rstrip('\n') for line in open(geneID)]
	sys.stdout = open(out_gff, 'w')
	for ID in geneIDs:
		children = db.children(ID, featuretype='exon')
		for feature in children:
			print(feature, "geneID=" + ID , "parenttype=" + ft , sep=";")
	sys.stdout.close()


def Gff2BedMod(bedfile):
	'''This function modifies the gff2bed output by dropping uncessary fields and 
	rearranging the bedfile such that it can be bedtools merged.

	Note, original gff2bed output raises error when bedtoold merged.  

	e.g
	gff2bed output: 
	2L      8192    8589    eID=FBgn0031208 .       +       FlyBase exon    .       Parent=FBtr0300690;geneID=FBgn0031208;parenttype=mRNA

	Gff2BedMod mod output:
	2L      8192    8589    Parent=FBtr0300690;geneID=FBgn0031208;parenttype=mRNA   .       +
	'''

	bed = pd.read_csv(bedfile, sep='\t', header=None, index_col=None)
	bed[3] = bed[9]
	bed.drop([6,7,8,9], axis =1, inplace = True)
	bed.to_csv(bedfile, sep='\t', header=None, index=False)



def CisNatBedMod(bedfile):
	
	'''This function converts the output of 'bedtools intersect -wb' 
	to a .bed format that has information of both strands in [3].
	Intersected strand details separated by underscore. 
	e.g plus_minus, mRNA_other
	'''
	
	cisnat_info = pd.read_csv(bedfile, sep='\t', header=None, index_col=None)
	cisnat_info[3] = cisnat_info[3].astype(str) + "_" + cisnat_info[9].astype(str)
	cisnat_info[5] = "." 
	cisnat_info.drop([6,7,8,9,10,11], axis =1, inplace = True)
	cisnat_info.to_csv(bedfile, sep='\t', header=None, index=False)


def RmBedMod(in_rmfile, out_bedfile, rm_ft):
	
	'''This function takes in a .RMout file and extracts all features 
	of each rm_ft and converts them to .bed format
	'''

	rm_file = pd.read_csv(in_rmfile, sep = '\s+', header = None, index_col= None, skiprows = 2)
	rm_bed  = rm_file[[4,5,6,10,7,8]]
	rm_bed[8] = rm_bed[8].str.replace('C', '-')
	rm_bed[7] = '.'
	rm_feature = rm_bed.loc[rm_bed[10].str.contains(rm_ft)]
	rm_feature.to_csv(out_bedfile, sep='\t', header=None, index=False)


def MakeGenomeMapLog(mapfile,unmapfile,lib,outlog):

	'''This function takes in both the mapped.txt and umapped.fasta files from genome mapping 
	and generates a maplog file
	'''
	 
	mapped = pd.read_csv(mapfile,sep='\t', header=None, index_col=None)
	unmapped = pd.read_csv(unmapfile,sep='\t', header=None, index_col=None)
	mapped['is_multiple_mapper'] = mapped[0].duplicated(keep = False)
	single_map_ID =mapped.loc[mapped['is_multiple_mapper'] == False][0].str.split('_', expand = True)[2].astype('int')
	multi_map_ID = mapped.loc[mapped['is_multiple_mapper'] == True][0].drop_duplicates().str.split('_', expand = True)[2].astype('int')
	uniq_map_ID = mapped[0].drop_duplicates().str.split('_', expand = True)[2].astype('int')
	unmap_ID = unmapped.loc[unmapped[0].str.startswith('>')][0].str.split('_', expand = True)[2].astype('int') #converting to series.

	table = {'Type' : ['Reads', 'Read_Counts'],
			'Mapped' : [len(uniq_map_ID),uniq_map_ID.sum()],
			'Unmapped' : [len(unmap_ID), unmap_ID.sum()],
			'PercentMapped': ['{}'.format(round(float(len(uniq_map_ID)) / float(len(uniq_map_ID) + len(unmap_ID)),3) * 100),'{}'.format(round(float(uniq_map_ID.sum()) / float(uniq_map_ID.sum() + unmap_ID.sum() ),3) * 100)],
			'UniqueMappers': [len(single_map_ID),single_map_ID.sum()],
			'MultipleMappers' : [len(multi_map_ID), multi_map_ID.sum()],
			'PercentMultipleMappers': ['{}%'.format(round(float(len(multi_map_ID)) / float(len(single_map_ID) + len(multi_map_ID)),3) * 100),'{}%'.format(round(float(multi_map_ID.sum()) / float(single_map_ID.sum() + multi_map_ID.sum()),3) * 100)],
			}

	out_df = pd.DataFrame(table, columns=['Type','Mapped','Unmapped','PercentMapped','UniqueMappers','MultipleMappers','PercentMultipleMappers'])
	print(lib, ' mapping statistics:')
	print(out_df)
	out_df.to_csv(outlog, sep='\t')


def SortMapFiles(filename):
	'''This function sorts the names of files in os.listdir()
	when gathering counts for each library. It ensures that
	mapped files with 'Index{num}' are sorted based on {num} 
	'''
	tag = filename.split('_')[1]
	if tag.startswith('index'):
		return int(tag.split('x')[1])
	else:
		return 0

def GatherBowtieCounts(lib_list, index,bowtie_out_dir, counts_out_dir):

	''' This function outputs the total number of reads mapped 
	to a given bowtie index by extracting counts from a bowtie output file.
	'''

	holder_list = []
	for lib in lib_list:
		file_path = pjoin(bowtie_out_dir, lib + '_' + index + '_mapped.txt')
		if os.path.getsize(file_path) > 0:
			mapped =  pd.read_csv(file_path,sep='\t', header=None, index_col=None)
			mapped_counts = mapped[0].str.split('_', expand = True)[2].astype('float')
			mapped_counts_sum = mapped_counts.sum()
		else:
			mapped_counts_sum = 0
		holder_list.append([lib,mapped_counts_sum])
	df = pd.DataFrame(holder_list,columns= ['Libraries', index])
	df.to_csv(pjoin(counts_out_dir, index + '_rawcounts.txt'), sep='\t', index=False)

def GatherLengthfilterBowtieCounts(lib_list, index,bowtie_out_dir, counts_out_dir, min_value, max_value):

	''' This function outputs the total number of reads that are less/greater than equal a max/min length
	mapped  to a given bowtie index by extracting counts from a bowtie output file.
	'''
	
	holder_list = []
	for lib in lib_list:
		file_path = pjoin(bowtie_out_dir, lib + '_' + index + '_mapped.txt')
		if os.path.getsize(file_path) > 0:
			mapped =  pd.read_csv(file_path,sep='\t', header=None, index_col=None)
			##mapped = mapped.loc[mapped[0].str.contains('21nt')]
			mapped.length = mapped[0].str.split('_', expand = True)[3].astype('str')
			mapped.length = pd.to_numeric(mapped.length.str.split('nt', expand=True)[0])
			mapped= mapped[(mapped.length >= min_value) & (mapped.length <= max_value)]
			if mapped.empty:
				mapped_counts_sum = 0
			else:
				mapped_counts = mapped[0].str.split('_', expand = True)[2].astype('float')
				mapped_counts_sum = mapped_counts.sum()
		else:
			mapped_counts_sum = 0
		holder_list.append([lib,mapped_counts_sum])
	df = pd.DataFrame(holder_list,columns= ['Libraries', index])
	df.to_csv(pjoin(counts_out_dir, index + '_rawcounts.txt'), sep='\t', index=False)