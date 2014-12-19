### 1\_data\_preproc.R

Pre-processing of 1000 Genomes and Gourraud et al (2014) datasets for later comparison of genotypes and allele frequencies

### 1\_data\_preproc\_figs.R

Generates figures with data from 1\_data\_preproc.R

### 2\_genotype\_mismatch.R

Compare genotypes in ARS exons of HLA between 1000 Genomes and Gourraud et al (2014)

### 2\_genotype\_mismatch\_figs.R

Generates figures with data from 2\_genotype\_mismatch.R

### 3\_af\_accuracy.R

Compare SNP allele frequencies in 1000 Genomes and Gourraud et al (2014) dataset

### 3\_af\_accuracy\_figs.R

Generates figures with data from 3\_af\_accuracy.R

### 4\_coverage.R

Relate genotype mismatch and frequency deviation between 1000G and PAG2014 to coverage

### 4\_coverage\_figs.R

Generates figures with data from 4\_coverage.R

### get\_coverage.py

Runs samtools to get coverage from 1000G bam files

### hla\_diff\_ref.R 

Output is a distance matrix of all HLA Class I alleles from the allele of the reference genome in a given window

### hla\_ref.R

Creates list of sequences of exons 2 and 3 (only exon 2 for class II) of the HLA alleles present in the reference human genome

### make\_consensus.R 

Creates consensus sequences at the 4-digit resolution for all HLA Class I and II alleles using IMGT alignment data v3170

### pag\_consensus.R 

Generates sequences for HLA alleles in Gourraud et al (2014) dataset using IMGT alignment data v2260

### pg2014\_vs\_erlich.R

Compare HLA alleles in PG2014 and Erlich goldstandard dataset
