#!usr/bin/python

###################
# Get coverage from bam file:
# chr1 pos1 chr2 pos2
###################

# command: python 01.py path/to/bam/file.bam
# input: path to bamfile (bai file should be in same folder)

import pysam
import re
import os
import sys
import numpy as np

if len(sys.argv) < 2:
	print '!!! ERROR: Please give the path to the bamfile as input !!!'
	os._exit(0)

bamfile_path = sys.argv[1]
os.mkdir('CXC_output')
header = '''samtools view -H ''' + bamfile_path + ''' > CXC_output/discordant.sam'''
os.system(header)
cmd_dis = '''samtools view ''' + bamfile_path + ''' | awk '$7 != "="' >> CXC_output/discordant.sam'''
os.system(cmd_dis)
sam_to_bam = '''samtools view -bS CXC_output/discordant.sam > CXC_output/discordant.bam'''
os.system(sam_to_bam)
os.system('samtools index CXC_output/discordant.bam')

#the discordant read bam file is read in 
bamfile = pysam.AlignmentFile("CXC_output/discordant.bam", "rb")

chromosomes = ['1','2','3','4','5','6','7','8','9','10','11', '12', '13', '14', '15', 
				'16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']

rep = {'X':'23', 'Y':'24', 'M':'25'}
rep = dict(rep.iteritems())
pattern = re.compile("|".join(rep.keys()))

for i in chromosomes: #for each chromosome the coverage information is saved in its own file
	chr = "chr" + i
	filename = 'CXC_output/'+ chr + '.txt'
	f = open(filename,'w')

	for column in bamfile.pileup(chr): #Loops through every position of the current chromosome
		#then get the information for all the reads at this position and their mate reads
		cov = bamfile.count_coverage(chr, column.pos, column.pos)
		print cov 
		cov_ar = np.zeros((cov,4))
		read_nr = 0

		i = pattern.sub(lambda m:rep[re.escape(m.group(0))],i)
		for read in column.pileups:
			cov_ar[read_nr,0] = i
			cov_ar[read_nr,1] = column.pos
			cov_ar[read_nr,2] = pattern.sub(lambda m: rep[re.escape(m.group(0))],str(list(read.alignment.next_reference_name)[3]))
			cov_ar[read_nr,3] = read.alignment.next_reference_start
			read_nr += 1

		chr_ar = np.bincount(cov_ar)
		if np.max(chr_ar)/cov >= 0.9:
			real_reads = cov_ar[cov_ar[:,2] == np.argmax(chr_ar)]
			outfile.write(real_reads)

	f.close()
