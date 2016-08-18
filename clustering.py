#!usr/bin/python

import numpy as np
import pandas as pd


#function that takes in a pd-series and max difference between two positions
def get_regions(series, windowsize):
	#find the start of a new region
	start_points = series[series.diff() > windowsize]
	start = np.array(start_points.index)
	start = np.append(start, len(series))
	return np.insert(start,0,0) #return np array with the beginning points of all regions in this series 


#open all the filtered chromosome files with format chr1 pos1 chr2 pos2 (one line per read)
chromosomes = ['1','2','3','4','5','6','7','8','9','10','11', '12', '13', '14', '15', 
				'16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']
outfile = open("outfile04.txt","w")


for chromosome in chromosomes:

	#open up the current file as a pandas dataframe
	infile = 'chr' + chromosome + '_3.txt'
	file_df = pd.read_table(infile, names = ["chr1", "pos1", "chr2", "pos2"], header= None)
	
	#extract the column containing position on the first chromsome
	pos1 = file_df["pos1"]
	#retrieve all the start points of new regions for this chromosome
	pos1_list = get_regions(pos1, 50)

	#loop through every region 
	for i in range(len(pos1_list)-1):
		#extract all columns within this region and find start and end of the region
		region_df = file_df[pos1_list[i]:pos1_list[i+1]]
		region_df = region_df.reset_index(drop=True)
		min_pos1 = region_df["pos1"][0]
		max_pos1 = region_df["pos1"][len(region_df.index)-1]
		region_len_1 = max_pos1 - min_pos1 
		coverage_pos1 = len(region_df.index)
		#check whether the region is at least 10bp long (1 read is 150bp) !!! ADD FILTER FOR MINIMUM COVERAGE
		if (region_len_1 > 10) and (coverage_pos1 > 100):
			#get all the values for the second position and sort them (reindex the sorted series)
			pos2 = region_df["pos2"].sort_values()
			pos2 = pos2.reset_index(drop = True)

			#also find the startpoints of new regions if available
			pos2_list = get_regions(pos2, 450)

			#loop through all the startpoints for this region and extract start and endpoints for this region
			for n in range(len(pos2_list)-1):
				min_pos2 = pos2[pos2_list[n]] 
				max_pos2 = pos2[pos2_list[n+1]-1]
				region_len_2 = max_pos2 - min_pos2 
				coverage_pos2 = len(pos2.index)

				#check whether region is 10bp long !!! ADD FILTER FOR MINIMUM COVERAGE
				if (region_len_2 > 10) and (coverage_pos2 >= 20):

					outfile.write('\t'.join(np.array([region_df["chr1"][0] , min_pos1, max_pos1, region_df["chr2"][0], min_pos2, max_pos2, coverage_pos2, '\n'], dtype = str)))






