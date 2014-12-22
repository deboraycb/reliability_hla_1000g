# pag_consensus.R

# consensus sequences for Gourraud et al (2014) data HLA alleles

# data directory
data_dir <- "../data/"

# download IMGT alignment data version 2.26.0
alig_dir <- paste0(data_dir, "Alignments_Rel_2.26.0/")
if(!file.exists(alig_dir)) {
    download.file(paste0("ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/",
                         "previous_releases/prev_nomenclature/Alignments_Rel_2.26.0.zip"),
                  destfile = "alignments.zip")
    unzip("alignments.zip", exdir = alig_dir)
    file.remove("alignments.zip")
}

# PAG2014 data
pag2014 <- read.table(paste0(data_dir, "mhc.tab"), sep = "\t", 
		   header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# Discard CHD individuals, which are not in 1000G Phase I
pag2014 <- pag2014[- grep("Chinese from Denver-Colorado, USA$", pag2014$Population), ]

# Allele names history file
if(! file.exists(paste0(data_dir, "Allelelist_history.txt"))) {
  download.file("ftp://ftp.eimb.relarn.ru/hla/Allelelist_history.txt",
		destfile = paste0(data_dir, "Allelelist_history.txt"))
}

# read allele names history table
allele_hist <- read.table(paste0(data_dir, "Allelelist_history.txt"),
			 sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# load function to parse IMGT alignment
source("hla_processnuc.R")

# function that takes a hla allele name in old nomenclature
## and returns the name as in a new version of the nomenclature
hla_subnames <- function(x, new.ver = 3170, old.ver = 2260, 
			 allele.hist = allele_hist)
{
  newnames <- allele.hist[, paste0("X", new.ver)]
  oldnames <- allele.hist[, paste0("X", old.ver)]
  newnames[which(oldnames == x)]
}

for(loc in c("A", "B", "C", "DRB1", "DQB1")) {
    
    # erase the digit at the end of DRB1 and DQB1 names
    ## IMGT file names don't have the digit 1
    locus <- gsub("[0-9]$", "", loc)

    # read alignment data
    res <- hla_processnuc(locus, alig.dir = alig_dir)
    ref_al <- res[[1]]
    dat_df <- res[[2]]

    # replace names with newer nomenclature
    ref_al <- hla_subnames(ref_al)
    dat_df[, 1] <- sapply(1:nrow(dat_df), function(i) hla_subnames(dat_df[i, 1])) 

    # make a matrix with all allele sequences at the desired exons
    # If class I, extract exons 2 and 3
    if(loc %in% c("A", "B", "C")) ex.end = 3
    # If class II, extract exon 2
    if(loc %in% c("DRB1", "DQB1")) ex.end = 2
    ##all_seqs <- hla_matrix(dat_df, ex.start = 2, ex.end = 3)
    #create list "nseq", where each entry will be an allele
    nseq <- list()
    for(i in 1:nrow(dat_df)) {
    nseq[[i]] <- strsplit(dat_df[i,2], "\\|")[[1]]
    }
    nseq <- lapply(nseq, "[", 2:ex.end)

    # transform each element in the list into a char vector 
    string_seq <- lapply(lapply(nseq, strsplit, ""), unlist)

    # transform the list of matrices into a matrix
    ##in which each row will be an allele.
    # Number of positions in the sequence
    bp <- length(string_seq[[1]])    
    all_seqs <- matrix(NA, nrow = length(nseq), ncol = bp)
    rownames(all_seqs) <- dat_df$alleles
    for(i in 1:length(string_seq)) all_seqs[i,] <- string_seq[[i]]


    # for positions with "-", take the nucleotide in the reference allele
##all_seqs <- hla_undash(all_seqs, ref_al)
    ref <- which(rownames(all_seqs) == ref_al)

    for(i in 1:ncol(all_seqs)) {
        all_seqs[, i][all_seqs[, i] == "-"] <- all_seqs[ref, i]
    }

    # Trim all sequences to the position which are not "."
    ##in the reference allele.
##all_seqs <- hla_rmdels(all_seqs)
    all_seqs <- all_seqs[, all_seqs[ref, ] != "."]
    colnames(all_seqs) <- 1:ncol(all_seqs)

    # extract all allele names for this locus
    alleles <- character()
    pos <- grep(paste0("\\.", loc, "\\."), names(pag2014))
    for(i in pos) alleles <- c(alleles, pag2014[, i])
    alleles <- unique(alleles)

    # erase "XX" in the allele name "14XX"
    ##sequence for this allele will be a consensus of all alleles 14
    alleles <- gsub("X", "", alleles)

    # create a pattern to be used in regex to match allele names 
    patt <- sapply(1:length(alleles), function(i) 
         strsplit(alleles[i], "/"))

    prefix <- paste0("^", loc, "\\*")
    for(i in 1:length(patt)) {
    patt[[i]] <- paste0(prefix, patt[[i]], "(:|[NQLS]$|$)", collapse = "|")
    }

    patt <- unlist(patt)

    # create a matrix in which we will place the consensus sequences
##consens <- hla_consensus(all_seqs, alleles, pattern = patt)
    consens <- matrix(NA, nrow = length(alleles), ncol = ncol(all_seqs))
    rownames(consens) <- alleles
    colnames(consens) <- 1:ncol(all_seqs)
  
    for(k in 1:length(patt)) {
        # list of all alleles in the 4-digit group
        tmp <- all_seqs[grep(patt[k], rownames(all_seqs)), , drop = FALSE]
        # take the different nucleotides found in each position
        umb <- apply(tmp, 2, unique)
        # for each position, consider only nucleotide symbols,
        ##"*" for "not sequenced" and "." for deletion
        for(i in 1:length(umb)) {
        umb[[i]] <- umb[[i]][grep("[AGCT\\*\\.]", umb[[i]])]
    }
    # collapse each umbiguous position with "/"
    conseq <- sapply(umb, paste, collapse = "/")
    # put this sequence in the output matrix consens
    consens[k, ] <- conseq
    }

# assign a name for the output matrix paste the locus name
assign(paste0("cons_seq_pag_", loc), consens)
}

# save a list of output matrices
save(list = ls()[grep("cons_seq_pag_", ls())],
 file = paste0(data_dir,"cons_seq_pag.rda"))
