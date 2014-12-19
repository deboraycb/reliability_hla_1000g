################################################################################
#
# 5_refbias.R
#
# Direction of frequency deviation & reference allele bias
#
################################################################################

# 0. Load previous data

load("4_coverage.rda")

#-------------------------------------------------------------------------------
# 1. Quantify deviations in direction of ref or alt allele

# Count number of populations with higher, lower or equal frequency of
# reference allele in 1000G relative to PAG2014 in each SNP
higherfreq <- list()
lowerfreq <- list()
equalfreq <- list()

# threshold to define "equal" frequency of reference allele
wellthr <- 0.01

for(l in loci){
    higherfreq[[l]] <- apply(freqdev[[l]][-1,] > freqthr, 2, sum)
    lowerfreq[[l]] <- apply(freqdev[[l]][-1,] < (-1 * freqthr), 2, sum)
    equalfreq[[l]] <- apply(abs(freqdev[[l]][-1,]) < wellthr, 2, sum)
}

# Define a threshold for the number of populations with over or underestimation
# of reference allele to call a given site over or underestimated
npopthr <- 2
overest <- lapply(higherfreq, function(x) names(x[x>npopthr]))
underest <- lapply(lowerfreq, function(x) names(x[x>npopthr]))
wellest <- lapply(equalfreq, function(x) names(x[x>npopthr]))

sapply(overest,length)
sapply(underest,length)
sapply(wellest,length)

#-------------------------------------------------------------------------------
# 2. List of sites with poorly estimated frequencies

badfreq <- list()
for(l in loci){
    badfreq[[l]] <- as.numeric(c(overest[[l]],underest[[l]]))
}

# List of hg19 coordinates of sites with poorly estimated frequencies
for (l in loci){
    write.table(exars_hg19[[l]][badfreq[[l]]], 
                file = paste0(data_dir, l, "_untrusted.txt"),
                quote = F, row.names = F, col.names = F)
}

#-------------------------------------------------------------------------------
# 3. Get difference from reference allele for windows around each site

# Function to call script to find difference from ref allele at a window
# centered on a snp
hla_diff_ref <-function(loc,site,win)
{
    # create the command string and call the command using system()
    command=paste("Rscript hla_diff_ref_pierre.R", loc, site, win)
    cat(command,"\n")
    try(system(command))
    # load workspace containing matrix of differences to ref allele
    return(paste0("matdist",loc,site,".RDa"))
}

# matrices of distance to reference HLA allele will be stored in list with one
# element per HLA gene containing one element per SNP-centered window
w=50
hladist <- list()
for (l in loci){
    hladist[[l]] <- matrix(nrow = 2*nrow(pag2014[[l]]),
                           ncol = length(polim.kg[[l]]),
                           dimnames = list(NULL, as.character(polim.kg[[l]])))
    for (s in polim.kg[[l]]){
        load(hla_diff_ref(l,as.numeric(s),w))
        hladist[[l]][,as.character(s)] <- matdist
    }
}



#----------------------------------------------------------------------------------
# 4. Create boolean matrix to classify each SNP-allele as REF or !REF

isref.fun <- function(x,r){
    res <- matrix(nrow=nrow(x),ncol=ncol(x),dimnames=dimnames(x))
    for (s in 1:ncol(x)){
        for(i in 1:nrow(x)){
            res[i,s] <- grepl(r[s], strsplit(x[i,s], "/"))
        }
    }
    return(res)
}

isref <- list()
for(l in loci){
    all_seqs_pag <- rbind(pag2014[[l]][, polim.kg[[l]], 1],
                          pag2014[[l]][, polim.kg[[l]], 2])
    rs <- refseq[[l]][polim.kg[[l]]]
    # substitute non genotyped regions by refseq
    all_seqs_pag <- t(apply(all_seqs_pag,1,function(x){
                            nonseq <- grep("^\\*$", x)
                            if(length(nonseq)>0) x[nonseq] <- rs[nonseq]
                            return(x)
    }))
    isref[[l]] <- isref.fun(all_seqs_pag,  rs)
}

#----------------------------------------------------------------------------------
# 5. Separate groups of over, under and well estimated sites and
# groups of reference and alternative allele bearing windows

refdist_overest <- list()
altdist_overest <- list()
for (l in loci[c(-3,-5)]){
    for (s in overest[[l]]){
        if((as.numeric(s) > 25 & as.numeric(s) < 246) |
           (as.numeric(s) > 295 & as.numeric(s) < 520)){
            refdist_overest[[l]] <- c(refdist_overest[[l]],
                                      hladist[[l]][,s][isref[[l]][,s]])
            altdist_overest[[l]] <- c(altdist_overest[[l]],
                                      hladist[[l]][,s][!isref[[l]][,s]]-1)
        }
    }
}

refdist_underest <- list()
altdist_underest <- list()
for (l in loci[-5]){
    for (s in underest[[l]]){
        if((as.numeric(s) > 25 & as.numeric(s) < 246) |
           (as.numeric(s) > 295 & as.numeric(s) < 520)){
            refdist_underest[[l]] <- c(refdist_underest[[l]],
                                       hladist[[l]][,s][isref[[l]][,s]])
            altdist_underest[[l]] <- c(altdist_underest[[l]],
                                       hladist[[l]][,s][!isref[[l]][,s]]-1)
        }
    }
}

refdist_wellest <- list()
altdist_wellest <- list()
for (l in loci[-5]){
    for (s in wellest[[l]]){
        if((as.numeric(s) > 25 & as.numeric(s) < 246) |
           (as.numeric(s) > 295 & as.numeric(s) < 520)){
            refdist_wellest[[l]] <- c(refdist_wellest[[l]],
                                      hladist[[l]][,s][isref[[l]][,s]])
            altdist_wellest[[l]] <- c(altdist_wellest[[l]],
                                      hladist[[l]][,s][!isref[[l]][,s]]-1)
        }
    }
}



mean(unlist(refdist_overest))
mean(unlist(altdist_overest))
wilcox.test(unlist(refdist_overest),unlist(altdist_overest),alternative='l')

mean(unlist(refdist_wellest))
mean(unlist(altdist_wellest))
wilcox.test(unlist(refdist_wellest),unlist(altdist_wellest),alternative='l')

mean(unlist(refdist_underest))
mean(unlist(altdist_underest))
wilcox.test(unlist(refdist_underest),unlist(altdist_underest),alternative='l')

##############################################################################################
save.image("5_refmapbias.rda")
