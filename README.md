# CXC_new
This is a collection of scripts to find interchromosomal translocations from a bamfile. 
It will shortly be updated, this is provisional code. 

So far the code consists of 5 scripts:

(1) coverage.py:
    - Creates a bamfile consisting only of discordant reads from the original bam file
    - Uses the discordant bam file to access the coverage information at every position in the genome
    - checks whether 90% of discordant reads are mapped to the same chromosome and prints them out
    
(2) clustering.py:
    - finds regions of translocations and sums them up
  
(3) reverse_check.py:
    - checks, whether a translocation region was identified as such independently on both chromosomes involved
    
(4) exact_bp.py:
    - finds the exact regions for every breakpoint by extracting soft-clipped reads from the original bamfile
    
--> At the moment a wrapper function is being created to simplify the process of using this script.
