#!usr/bin/python

import pandas as pd
import numpy as np

outfile = open('output_final.txt', 'w')
full_df = pd.read_table('outfile_04.txt', header = None, names = ['chr1','min1','max1','chr2','min2','max2','cov','x'])
full_df = full_df.drop('x',1) 
full_df = full_df[full_df['chr1'] != full_df['chr2']]

reg_1 = full_df[['chr1', 'min1', 'max1']]
reg_2 = full_df[['chr2', 'min2', 'max2']]

for index, row in reg_1.iterrows():
	reg_2_chr = reg_2[reg_2['chr2'] == row['chr1']]
	reg_2_chr = reg_2_chr.sort_values('min2')
	min_fit = reg_2_chr['min2'] == row['min1'] 
	#max_fit = reg_2_chr['max2'] <= row['max1'] 
	#fit_rows = reg_2_chr[min_fit & max_fit]
	fit_rows = reg_2_chr[min_fit]
	i_of_matches = fit_rows.index
	if index in full_df.index:
		ori_row = full_df.ix[index]
		for i in i_of_matches:
			if i in full_df.index:
				partner_row = full_df.ix[i]
				if (ori_row['chr2']==partner_row['chr1']) and (ori_row['min2']== partner_row['min1']):
				#if (ori_row['chr2']==partner_row['chr1']) and (ori_row['min2']>= partner_row['min1']) and (ori_row['max2']<= partner_row['max1']):
					string1 = ''
					string2 = ''
					for word in ori_row:
						string1 = string1 + str(word) + '\t'
					for word in partner_row:
						string2 = string2 + str(word) + '\t'
					cov = str(ori_row['cov']+ partner_row['cov'])
					new_line = string1 + string2 + cov + '\n'
					outfile.write(new_line)
					full_df = full_df.drop(i, 0)
