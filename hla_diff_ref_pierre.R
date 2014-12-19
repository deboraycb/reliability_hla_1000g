#hla_diff_ref.R

# distance matrices of hla class I allele groups against the allele of the
# reference human genome

#load object with reference sequence at ARS exons
load("hla_ref.RData")

#load consensus sequences of Pierre dataset
source("1_data_preproc.R")

# get locus, SNP position and window size from the command line arguments
locus = as.character(commandArgs(trailingOnly=TRUE)[1])
snp = as.numeric(commandArgs(trailingOnly=TRUE)[2])
window_size = as.numeric(commandArgs(trailingOnly=TRUE)[3])

# define allele of the reference genome and get its sequence
refseq = refseq[[locus]]
all_seqs = rbind(pag2014[[locus]][,,1], pag2014[[locus]][,,2])

# define window
window_range = (snp-(window_size %/% 2)):(snp+(window_size %/% 2))

# if window exceeds the length of the sequence, stop and launch an error msg
if(any(!window_range%in%1:ncol(all_seqs))) {
##    stop(message("ERROR: chosen window exceeds sequence length, which is ",
##		 ncol(all_seqs), "."))
    window_range <- window_range[window_range%in%1:ncol(all_seqs)]
    cat(window_range,'\n')
}

# trim the sequences to the window size
all_seqs_wind = all_seqs[,window_range]

#all_seqs_wind = unique(all_seqs_wind)
refseq_wind = refseq[,window_range,drop=FALSE]

# substitute non genotyped regions by refseq
all_seqs_wind <- t(apply(all_seqs_wind,1,function(x){
      nonseq <- grep("^\\*$", x)
      if(length(nonseq)>0) x[nonseq] <- refseq_wind[,nonseq]
      return(x)
}))


# create a distance matrix of each allele and the reference allele
matdist = t(apply(all_seqs_wind, 1, function(x) {
                 res <- numeric(ncol(refseq_wind))
                    for (i in 1:ncol(refseq_wind)){
                        res[i] <- mean(sapply(refseq_wind[,i], grepl,
                                              strsplit(x[i], '/')[[1]]))
                    }
                 return(sum(1 - res))
}))

# save file and finish program
outfile_name = paste0("matdist", locus, snp, ".RDa")
outfile = paste0(getwd(), "/", outfile_name)
save(matdist, file=outfile)
message("DONE! Output file ", outfile_name,
        " was saved to your current directory")
