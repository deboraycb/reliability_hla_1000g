################################################################################
#
# 3_af_accuracy.R
#
# Allele frequency deviation & reference allele bias
#
################################################################################

# source script with PAG2014 data, 1000G data and genotype comparisons
source("2_genotype_mismatch.R")

#-------------------------------------------------------------------------------
# 1. Frequency of REF allele in 1000G and in PAG2014

# source script to get sequences of HLA alleles of the reference human genome
source("hla_ref.R")

# Function to count number of ref alleles in 1000G and PAG2014 data
countref <- function(kgl, pagl, refl, pop){
    # rbinds the 2 alleles at each locus generating an array of
    # (2alleles*n_pop_individuals) lines by (intersect_loci) columns
    kg = rbind(kgl[pop,,1], kgl[pop,,2])
    pag = rbind(pagl[pop,,1], pagl[pop,,2])
    cnt = matrix(dimnames = list(c("1000G", "PAG2014"), colnames(kg)),
                 nrow = 2, ncol = ncol(kg))
    for (i in 1:ncol(kg)){
        cnt[1,i] <- sum(kg[,i] == refl[i])
        cnt[2,i] <- sum(pag[,i] == refl[i])
    }
    return(cnt)
}

# List of individuals in each population
pop <- list()
for (p in unique(pg$Population)){
    pshort <- substr(p, 16, 18)
    pop[[pshort]] <- pg$Subject[pg$Population == p]
}
pop <- c(list(global = pg$Subject), pop)

# lists of counts and frequency of reference allele per locus per population
nref <- list()
reffreq <- list()
for (l in loci){
    nref[[l]] <- array(dim = c(length(pop), ncol(kg[[l]]), 2),
                       dimnames = list(names(pop), colnames(kg[[l]]),
                                       c("1000G","PAG2014")))
    reffreq[[l]] <- array(dim = c(length(pop), ncol(kg[[l]]), 2),
                          dimnames = list(names(pop), colnames(kg[[l]]),
                                          c("1000G","PAG2014")))
    for(p in names(pop)){
        nref[[l]][p,,] <- t(countref(kg[[l]], pag2014[[l]], refseq[[l]],
                                     pop[[p]]))
        reffreq[[l]][p,,] <- nref[[l]][p,,] / (2 * length(pop[[p]]))
    }
}

# Diference in frequency
freqdev <- mapply(function(l, s)
                  l[,s,"1000G"] - l[,s,"PAG2014"],
                  reffreq, polim.kg)


#-------------------------------------------------------------------------------
# 2. Deviation of frequency in 1000G from expected frequency in PAG2014

# root mean squared error and mean absolute error
rmse <- matrix(nrow = length(pop), ncol = length(loci),
                dimnames = list(names(pop), loci))
mae <- matrix(nrow = length(pop), ncol = length(loci),
                dimnames = list(names(pop), loci))
p='global'
freqthr <- 0.1

for(l in loci){
    for(p in names(pop)){
        # Calculate sum of squared differences
        ssres <- sum((reffreq[[l]][p, , "1000G"] - reffreq[[l]][p, , "PAG2014"]
                         ) ^ 2, na.rm = T)
        # Calculate root-mean-squared-error
        rmse[p,l] <- sqrt(ssres / length(polim.kg[[l]]))
        # Calculate absolute error
        ae <- abs(reffreq[[l]][p, , "1000G"] - reffreq[[l]][p, , "PAG2014"])
        mae[p,l] <- sum(ae, na.rm = T) / length(polim.kg[[l]])
    }
}

# RMSE over all genes and merging populations
reffreq.ovl <- data.frame(kg = unlist(sapply(reffreq,
                                             function(x) x[1,,"1000G"])),
                          pag2014 = unlist(sapply(reffreq,
                                                  function(x) x[1,,"PAG2014"])))

ssres.ovl <- sum((reffreq.ovl$kg - reffreq.ovl$pag2014) ^ 2, na.rm = T)
rmse.ovl <- sqrt(ssres.ovl / sum(sapply(polim.kg, length)))

# MAE over all genes and merging populations
mae.ovl <- sum(abs(reffreq.ovl$kg - reffreq.ovl$pag2014), na.rm = T) /
                            sum(sapply(polim.kg, length))
