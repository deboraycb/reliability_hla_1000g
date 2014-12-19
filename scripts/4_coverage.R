################################################################################
#
# 4_coverage.R
#
# Relate genotype mismatch and frequency deviation to coverage
#
################################################################################

# Source script with PG2014 data, 1000G data, genotype mismatch and frequency
# difference objetcs
source("3_af_accuracy.R")

#-------------------------------------------------------------------------------
# 1. Get coverage per site per individual

# create one file for each population only with individuals in both datasets
for (p in unique(indkg$POP)){
    write.table(indkg$IND[indkg$POP == p & indkg$IND %in% pier$Subject],
                quote = F, row.names = F, col.names = F,
                paste0(data_dir, "pops/overlap_kg_pag_", p, ".txt"))
}

# and run get_coverage.py to get coverage for exons 2 and 3 from A B C
# DRB1 and DQB1 for all individuals in Pierre's dataset from BAM files in the
# 1000G FTP
try(system("python ~/hla_tools/1kg/scripts/get_coverage.py"))

### Function that calls indersectbed from bedtools
BEDintersect <- function(functionstring = "intersectBed", bed1, bed2,
                         opt.string=""){
    #create temp files
    a.file = tempfile()
    b.file = tempfile()
    out = tempfile()
    options(scipen = 99) # not to use scientific notation when writing out
    #write bed formatted dataframes to tempfile
    write.table(bed1, file = a.file, quote = F, sep = "\t", col.names = F,
                row.names = F)
    write.table(bed2, file = b.file, quote = F, sep = "\t", col.names = F,
                row.names = F)
    # create the command string and call the command using system()
    command = paste(functionstring, "-a", a.file, "-b", b.file, opt.string,
                    ">", out, sep = " ")
    cat(command, "\n")
    try(system(command))
    res = try(read.table(out, header = F))
    if(class(res) == 'try-error'){
        res = NA
    }
    unlink(a.file); unlink(b.file); unlink(out)
    return(res)
}

bed <- list()
cover <- list()
for (l in loci){
    # create bed files for exons 2 and 3 of A B and C and exon 2 from DQB1 and
    # DRB1 with one position per line
    bed[[l]] <- cbind(6, exars_hg19[[l]][polim.kg[[l]]] - 1,
                      exars_hg19[[l]][polim.kg[[l]]])
    #create matrices to store coverage per site per individual
    cover[[l]] <- matrix(nrow = nrow(mism[[l]]), ncol = ncol(mism[[l]]),
                         dimnames = dimnames(mism[[l]]))
}

# intersect the bed files above with the bed files containing coverage in the
# 4th column, so we have coverage for each site in each line
for (i in rownames(mismatch)){
    cat(i,'\n')
    # BED file with coverage for the whole MHC region for individual i
    tempcovbed <- read.table(paste0(data_dir, "/coverage/coverage_", i, ".bg"))
    templ <- list()
    for (l in loci){
        # intersect bed files 1) with one position per line and 2) with coverage
        # in the 4th column
        templ[[l]] <- BEDintersect(bed1 = tempcovbed, bed2 = bed[[l]])
        # merge with original bed file so that positions that don't have
        # coverage values will be included (and 4th col will be NA)
        templ[[l]] <- merge(bed[[l]], templ[[l]], all.x = T)
        # substitute NA per 0, because positions that are not in the coverage
        # bed actually have 0 coverage (they were imputed)
        templ[[l]][is.na(templ[[l]])] <- 0
        # save only coverage vector in cover.a object
        cover[[l]][i,] <- templ[[l]][,4]
    }
}

#-------------------------------------------------------------------------------
# 2. Genotype mismatches vs. coverage

# 2.1. Function to test for difference in means of two groups, permuting the
# groups.
permtest <- function(x, group, n = 1000, alternative = ""){
    obsdif = mean(x[group]) - mean(x[!group])
    onetest <- function(x, group){
        newgroup <- sample(group)
        mean(x[newgroup]) - mean(x[!newgroup])
    }
    permdif = replicate(n, onetest(x, group))
    if(alternative == "greater"){
        cat("\n Test: mean in group is greater than in !group\n")
        p = mean(permdif > obsdif)
    }
    if(alternative == "less"){
        cat("\n Test: mean in group is less than in !group\n")
        p = mean(permdif < obsdif)
    }
    else{
        cat("\n Test: mean in group is different from !group\n")
        p = mean(abs(permdif) > abs(obsdif))
    }
    return(list(o = obsdif, e = permdif, p = p, n = n, alt = alternative))
}

# 2.2. Get list of sites discovered only on exome experiments (thus with high
# coverage) for exclusion

for (l in loci){
    command <- paste0("vcftools --gzvcf ", data_dir,
                      "ALL.chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz",
                      " --bed ", data_dir, tolower(l), "_ex2*.bed --recode",
                      " --recode-INFO SNPSOURCE --stdout |",
                      " grep -E 'SNPSOURCE\\=EXOME\\s' | cut -f 2 > ",
                      data_dir, "exome_positions_", tolower(l), ".txt")
    system(command)
}

# Transform list of site positions to ARS exons coordinates
exome_pos <- list()
for (l in loci){
    temp <- read.table(paste0(data_dir, "exome_positions_",
                                        tolower(l), ".txt"))
    exome_pos[[l]] <- data.frame(hg19 = temp,
                                 ars = which(exars_hg19[[l]] %in% temp$V1))
    exome_pos[[l]] <- cbind(exome_pos[[l]], polim = which(polim.kg[[l]] %in%
                                                  exome_pos[[l]][,"ars"]))
}

# coverage and mismatch only for LOWCOV or LOWCOV,EXOME SNPs
cover_lc <- list()
mism_lc <- list()
polim.kg_lc <- list()
for (l in loci){
    cover_lc[[l]] <- cover[[l]][, -exome_pos[[l]][,3] ]
    mism_lc[[l]] <- mism[[l]][, -exome_pos[[l]][,3] ]
    polim.kg_lc[[l]] <- polim.kg[[l]][ -exome_pos[[l]][,3] ]
}

# 2.3. Permutation to test if coverage in mismatched genotypes < coverage in 
# matched genotypes
covergt.perm <- permtest(x = unlist(cover_lc), group = unlist(mism_lc),
                         alternative = "less", n = 10000)

# 2.4. Mann-whitney to test difference in coverage
covergt.mw <- wilcox.test(unlist(cover_lc)[unlist(mism_lc)],
                          unlist(cover_lc)[!unlist(mism_lc)],
                          alternative = "l")

#--------------------------------------------------------------------------------
# 3. Departure from expected frequency vs. coverage

# 3.1. Measured by absolute difference in frequency per site, for lowcov sites

freqdev_lc <- mapply(function(l, s)
                  l[,s,"1000G"] - l[,s,"PAG2014"],
                  reffreq, polim.kg_lc)

# Mean coverage per site
meancov <- sapply(cover_lc, function(x) apply(x, 2, mean))

# Linear model
freqcov.lm <- lm(unlist(sapply(freqdev_lc, function(x) abs(x["global", ]))) ~ 
                 unlist(meancov))
freqcov.sum <- summary(freqcov.lm)

# 3.2. Categorize sites in i) well estimated; 2) poorly estimated frequencies
freqthr <- 0.1
poorfreq <- lapply(freqdev, function(x) abs(x["global", ]) >= freqthr)

coverfreq.perm <- permtest(x = unlist(meancov), group = unlist(poorfreq),
                           alternative="less", n = 10000)

################################################################################
save.image("4_coverage.rda")
