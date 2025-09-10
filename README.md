# Pabietivora
Contains the python script to conduct windows sliding analysis of minor allele frequency along the Phytophthora genome for production of Figure 1 in Feau et al. (unpublished)

Input is a vcf file of Illumina reads remapped on the draft genome of P. abietiivora (JBQXYJ000000000). The script will graph a distribution of minor allele frequencies (Figure 1A) and then generate a graph for each scaffold (>100 kbp) with i. average minor allele frequency; ii. number of SNP loci; and iii. average sequencing read in 1 kbp windows (Figure 1B).

WSalllelesFreq.py require the installation of samtools.

usage: python WSalllelesFreq.py
