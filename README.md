### 1\_data\_preproc.R

Pre-processing of 1000 Genomes and Gourraud et al (2014) datasets for later comparison of genotypes and allele frequencies

### 2\_genotype\_mismatch.R

Compare genotypes in ARS exons of HLA between 1000 Genomes and Gourraud et al (2014)

### 3\_af\_accuracy.R

Compare SNP allele frequencies in 1000 Genomes and Gourraud et al (2014) dataset

### 4\_coverage.R

Relate genotype mismatch and frequency deviation between 1000G and PAG2014 to coverage

### get\_coverage.py

Runs samtools to get coverage from 1000G bam files

### hla\_diff\_ref\_pierre.R 

Output is a distance matrix of all HLA alleles from the allele of the reference genome in a given window

### hla\_ref.R

Creates list of sequences of exons 2 and 3 (only exon 2 for class II) of the HLA alleles present in the reference human genome

### pag\_consensus.R 

Generates sequences for HLA alleles in Gourraud et al (2014) dataset using IMGT alignment data v2260
