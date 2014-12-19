# pag_consensus.R

# consensus sequences for Gourraud et al (2014) data HLA alleles

# loads functions in Rfunctionslib/
source("~/hla_tools/Rfunctionslib/hla_library.R")
hla_library()

# data directory
data_dir <- "~/hla_tools/1kg/data/"

# download IMGT alignment data
hla_download(alig.version = "2.26.0", out.dir = data_dir, prev.nom = TRUE)

# Pierre's data
pier <- read.table(paste0(data_dir, "mhc.tab"), sep = "\t", 
		   header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# Discard CHD individuals, which are not in 1000G Phase I
pier <- pier[- grep("Chinese from Denver-Colorado, USA$", pier$Population), ]

# Allele names history file
if(! file.exists(paste0(data_dir, "Allelelist_history.txt"))) {
  download.file("ftp://ftp.eimb.relarn.ru/hla/Allelelist_history.txt",
		destfile = paste0(data_dir, "Allelelist_history.txt"))
}

# read allele names history table
allele_hist <- read.table(paste0(data_dir, "Allelelist_history.txt"),
			 sep = "\t", header = TRUE, stringsAsFactors = FALSE)

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
  # If class II, extract exon 2
  if(loc %in% c("A", "B", "C")) {
    all_seqs <- hla_matrix(dat_df, ex.start = 2, ex.end = 3)
  }

  if(loc %in% c("DRB1", "DQB1")) {
    all_seqs <- hla_matrix(dat_df, ex.start = 2, ex.end = 2)
  }

  # for positions with "-", take the nucleotide in the reference allele
  all_seqs <- hla_undash(all_seqs, ref_al)

  # Trim all sequences to the position which are not "."
  ##in the reference allele.
  all_seqs <- hla_rmdels(all_seqs)
  colnames(all_seqs) <- 1:ncol(all_seqs)

  # extract all allele names for this locus
  alleles <- character()
  pos <- grep(paste0("\\.", loc, "\\."), names(pier))
  for(i in pos) alleles <- c(alleles, pier[, i])
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
  consens <- hla_consensus(all_seqs, alleles, pattern = patt)
  
  # assign a name for the output matrix paste the locus name
  assign(paste0("cons_seq_pag_", loc), consens)
}

# save a list of output matrices
save(list = ls()[grep("cons_seq_pag_", ls())],
     file = "~/hla_tools/1kg/data/cons_seq_pag.rda")
