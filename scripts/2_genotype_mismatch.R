################################################################################
#
# 2_genotype_mismatch.R
#
# Compare genotypes in 1000G and PAG2014
#
# Oct/2014
#
################################################################################
#----------------------------------------------------------
# Source script that preprocesses 1000G and PAG2014 datasets
source("1_data_preproc.R")

#----------------------------------------------------------
# 1. Comparing genotypes (using only polymorphic positions at 1000G)

# Function to compare genotypes at each position in one individual
mymatch <- function(kg, sanger){
    sanger.new <- matrix(nrow = nrow(sanger), ncol = ncol(sanger))
    pos <- vector("logical", nrow(kg))
    # for each position i
    for (i in 1:nrow(kg)){
        # grep KG allele or * in sanger allele;
        tab <- sapply(kg[i,], grepl, sanger[i,])
        # if allele 1 in KG is compatible to allele 1 in Sanger and allele 2 to
        # allele 2
        if (sum(tab[cbind(1:2, 1:2)]) == 2){
            # This is a match
            # Sanger alleles (which might be ambiguous) are replaced by KG
            # alleles in the sanger.new matrix
            sanger.new[i,] <- kg[i,]
            #and mism will be atributed F for that individual, in that position
            pos[i] <- F
        }
        # if allele 1 in KG is compatible to allele 2 in Sanger and allele 2 to
        # allele 1
        else if (sum(tab[cbind(1:2, 2:1)]) == 2){
            # This is also a match
            # Sanger alleles (which might be ambiguous) are replaced by KG
            # alleles in the sanger.new matrix
            sanger.new[i,] <- kg[i,]
            #and mism will be atributed F for that individual, in that position
            pos[i] <- F
        }
        else{
            #Otherwise, this is a mismatch;
            #sanger.new will get the original sanger allele
            sanger.new[i,] <- sanger[i,]
            #and mism will be attributed T
            pos[i] <- T
        }
    }
    return(list(pos = pos, sanger.new = sanger.new))
}


#Then apply function above to all individuals

#mism.[a|b|c|dqb|drb] will be logical matrices with individuals in rows and loci
#in columns; it will be T where there's a mismatch and F where there's a match.

#pier.new.[a|b|c|dqb|drb] will be a 'corrected' Pierre dataset. Ambiguities will
#be replaced by the allele at the KG, if there's one.

mism <- list()
pag2014.new <- list()
loci <- c("A", "B", "C", "DQB1", "DRB1")

for (l in loci){
    mism[[l]] <- matrix(nrow = nrow(kg[[l]]), ncol = length(polim.kg[[l]]),
                        dimnames = list(rownames(kg[[l]]),
                                        polim.kg[[l]]))
    pag2014.new[[l]] <- array(dim = c(nrow(kg[[l]]), length(polim.kg[[l]]), 2),
                           dimnames = list(rownames(pag2014[[l]]),
                                           polim.kg[[l]], c("a1","a2")))
    for (i in 1:nrow(kg[[l]])){
        temp <- mymatch(kg[[l]][i, polim.kg[[l]], ],
                        pag2014[[l]][i, polim.kg[[l]], ])
        mism[[l]][i,] <- temp$pos
        pag2014.new[[l]][i,,] <- temp$sanger.new
    }
}

mismatch <- do.call(cbind, mism)
polim.all.kg <- do.call(c, polim.kg)

#mismatch is concentrated in a few sites:
tot = 0
for (n in 1:ncol(mismatch)){
    tot = tot + sort(apply(mismatch, 2, sum), decreasing = T)[n]
    if (tot >= sum(mismatch) / 2){
        cat(tot, "out of", sum(mismatch), "mismatches (", tot / sum(mismatch),
            ") concentrated on", n, "out of", ncol(mismatch), "sites (",
            n / ncol(mismatch), ").\n")
        break
    }
    }
