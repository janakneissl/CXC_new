#!usr/bin/python

import pysam
import numpy as np
import re
import pandas as pd
import itertools

rep = {'23':'X', '24':'Y', '25':'M'}
rep = dict(rep.iteritems())
pattern = re.compile("|".join(rep.keys()))

#open the bamfile with original alignment and the tab-delimited file with the potential breakpoint regions
bamfile = pysam.AlignmentFile('transloci.bam', 'rb')
regions = pd.read_table('output_final_0908_3.txt', header = None, names = ['t1_rchr', 
	't1_rst', 't1_ren', 't1_mchr', 't1_mst', 't1_men', 't1_cov', 't2_rchr', 
	't2_rst', 't2_ren', 't2_mchr', 't2_mst', 't2_men', 't2_cov', 'tot_cov' ] )
outfile = open('110816_06_disc.txt', 'w')
header_list = ['1. Chromsome', '1. Breakpoint', '2. Chromosome', '2. Breakpoint', 'Total Reads', 'Discordant Reads', 'Softclipped Reads', '\n']
header = '\t'.join(header_list)
outfile.write(header)

#drop unimportant columns from the breakpoint regions file
regions = regions.drop(['t1_mchr','t1_mst', 't1_men', 't1_cov', 't2_mchr', 't2_mst', 't2_men', 'tot_cov'],1)

#for each potential breakpoint:
for index,region in regions.iterrows(): 
	#get the chromosomes between which the breakpoint is suspected
	fetch_chr_1 = 'chr' + str(region['t1_rchr'])
	fetch_chr_2 = 'chr' + str(region['t2_rchr'])
	fetch_chr_1 = pattern.sub(lambda m:rep[re.escape(m.group(0))],fetch_chr_1)
	fetch_chr_2 = pattern.sub(lambda m:rep[re.escape(m.group(0))],fetch_chr_2)

	#create a dataframe, wich will contain the reads in the two breakpoint regions
	df_region = pd.DataFrame(columns = ['ID', 'read_chr', 'read_start', 'CIGAR', 'mate_chr', 'mate_start'])

	#extract information from reads in bamfile for both regions
	for read in bamfile.fetch(reference=fetch_chr_1, start=region['t1_rst'],end= region['t1_ren']):
		d = {'ID':read.query_name,'read_chr':read.reference_name,'read_start':read.reference_start, 'CIGAR':read.cigarstring, 'mate_chr':read.next_reference_name, 'mate_start':read.next_reference_start}
		row = pd.DataFrame.from_records(d, index=[0])      
		df_region = df_region.append(row)
	for read in bamfile.fetch(fetch_chr_2, region['t2_rst'],region['t2_ren']):
		d = {'ID':read.query_name,'read_chr':read.reference_name,'read_start':read.reference_start, 'CIGAR':read.cigarstring, 'mate_chr':read.next_reference_name, 'mate_start':read.next_reference_start}
		row = pd.DataFrame.from_records(d, index=[0])      
		df_region = df_region.append(row)
	df_region = df_region.reset_index(drop=True)

	#each unique read will have a unique Sequence ID (this way mates are not counted double)
	IDs = df_region['ID'].unique()
	total_reads = len(IDs)
        
	#filter out the discordant reads and count the amount of unique discordant reads (sometimes, same read is mapped to multiple positions)
	discordant = df_region[df_region['mate_chr']!=df_region['read_chr']]
	discordant_IDs = len(discordant['ID'].unique())
	perc_discordant =float(discordant_IDs)/float(total_reads)
        
	#find the discordant reads that are discordant between the two chromosomes from the regions and count the number
	t1 = (discordant['mate_chr'] == fetch_chr_2) & (discordant['mate_start'] >= region['t2_rst']) & (discordant['mate_start'] <= region['t2_ren'])
	t2 = (discordant['mate_chr'] == fetch_chr_1) & (discordant['mate_start'] >= region['t1_rst']) & (discordant['mate_start'] <= region['t1_ren'])  
	correct_discordant = discordant[t1 | t2]
	correct_dis_IDs = len(correct_discordant['ID'].unique())
	perc_correct_dis = float(correct_dis_IDs)/float(total_reads)

	#extract the softclipped reads in the region to suggest an exact breakpoint (they have S in their CIGAR string)
	#only allow one S in cigar string to avoid reads soft-clipped on both ends
	#devide the soft clipped reads of both regions
	soft_clipped = correct_discordant[[(len(re.findall('S', x)) == 1) for x in correct_discordant['CIGAR']]]

	#for both regions a breakpoint is calculated by using the soft-clipped reads
	for chromosome in [fetch_chr_1, fetch_chr_2]:
		soft_clipped_r = soft_clipped[soft_clipped['read_chr']==chromosome]
		#at least two soft-clipped reads should lay in this region
		if len(soft_clipped_r.index) > 2: #differentiate between right and left clipped regions
			right_clipped_p = re.compile(r'[0-9]+S$')
			right_clipped = soft_clipped_r[[bool([re.search(right_clipped_p, x)]) for x in soft_clipped_r['CIGAR']]]
			right_add = np.array([int("".join(itertools.takewhile(str.isdigit, x))) for x in right_clipped['CIGAR']])
			right_breakpoints = right_add + np.array(right_clipped['read_start'])                       
			left_clipped_p = re.compile(r'^[0-9]+S')
			left_clipped = soft_clipped_r[[bool([re.search(left_clipped_p, x)]) for x in  soft_clipped_r['CIGAR']]]
			left_add = np.array([int("".join(itertools.takewhile(str.isdigit, x))) for x in left_clipped['CIGAR']])
			left_breakpoints = left_add + np.array(left_clipped['read_start'])

			if (len(right_breakpoints >=1)) | (len(left_breakpoints >=1)): 
				breakpoints = np.concatenate((right_breakpoints, left_breakpoints))
				breakpoint_dic = {}
				for key in np.unique(breakpoints):
					breakpoint_dic[key] = 0
				for breakpoint in breakpoints:
					breakpoint_dic[breakpoint] += 1
				breakpoint = max(breakpoint_dic, key = breakpoint_dic.get) #find the breakregion that occured most often 
				outstring = chromosome + '\t' + str(int(breakpoint)) + '\t'
				outfile.write(outstring)
		else: 
			outstring = chromosome + '\t' + 'None' + '\t'
			outfile.write(outstring)
	outstring_2= str(total_reads) + '\t' +str(discordant_IDs) +'\t' + str(len(soft_clipped))+'\n'
	outfile.write(outstring_2)
 





