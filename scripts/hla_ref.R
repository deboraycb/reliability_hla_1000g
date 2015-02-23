# hla_ref.R

# Get sequences of HLA alleles of the reference human genome

# data directory
data_dir <- "../data/"

# check if alignment data exists in directory
# if not, download it.
if (!file.exists(paste0(data_dir,"Alignments Rel_3160"))) {
    download.file("ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/previous_releases/Alignments_Rel_3160.zip",
		  destfile = paste0(data_dir, "alignments.zip"))
    unzip(paste0(data_dir, "alignments.zip"))
    file.remove(paste0(data_dir, "alignments.zip"))
}

refseq=list()

for (locus in c('A','B','C','DQB','DRB')){
    # read alignment data
    filename = paste0(getwd(), "/Alignments Rel_3160/", locus, "_nuc.txt")
    alig = scan(filename, what="character", sep="\n")
    if (locus == "DQB" | locus=="DRB") locus <- paste0(locus, "1")
    alig = alig[grep(paste0(locus, "\\*"), alig)]
    
    # remove blank spaces
    for(i in seq_along(alig)) {
        x = unlist(strsplit(alig[i], " "))
        alig[i] = paste(x[x!=""], collapse=" ")
    }
    
    # format data as data frame
    dat_df = data.frame(
        alleles=sapply(strsplit(alig, " "),"[",1),
        sequences=sapply(1:length(alig), function(i) paste0(
    	sapply(strsplit(alig[i]," "), "[", -1), collapse="")),
        stringsAsFactors=FALSE
    )
    
    # join the sequences of all exons of an allele in a single row 
    dat_df = aggregate(sequences ~ alleles, data=dat_df, paste0, collapse="")
    
    # extract exons 2 and 3
    # each allele will be an element of the list nseq
    nseq = list()
    if(grepl("D",locus)) {
        ex <- 2
    } else {ex <- 2:3}
    for(i in 1:nrow(dat_df)) nseq[[i]] = strsplit(dat_df[i,2], "\\|")[[1]][ex]
    names(nseq) = dat_df$alleles 
    
    # transform sequences into character vectors
    string_seq = lapply(lapply(nseq, strsplit, ""), unlist)
    
    # sequences' length
    bp = length(string_seq[[1]])    
    
    # format data as matrix
    all_seqs = matrix(NA, nrow=length(nseq), ncol=bp)
    rownames(all_seqs) = names(nseq)
    for(i in 1:length(string_seq)) all_seqs[i,] = string_seq[[i]]
    
    # get the nucleotide in the first allele for those positions with "-"
    ref = which(rownames(all_seqs)==unlist(strsplit(alig[1], " "))[1])
    for(i in 2:nrow(all_seqs)) {
        all_seqs[i,][all_seqs[i,]=="-"] = all_seqs[ref,][all_seqs[i,]=="-"]
    }
    
    # delete positions where all allele have a "."
    all_seqs = all_seqs[,all_seqs[ref,] != "."]
    colnames(all_seqs) = 1:ncol(all_seqs)
    
    # define allele of the reference genome and get its sequence
    if(locus=="A") ref="A*03:01:01:01"
    if(locus=="B") ref="B*07:02:01"
    if(locus=="C") ref="C*07:02:01:01"
    if(locus=="DQB1") ref="DQB1*06:02:01"
    if(locus=="DRB1") ref="DRB1*15:01:01:01"
    
    refseq[[locus]] = all_seqs[rownames(all_seqs)%in%ref,,drop=FALSE] 
    
}

save(refseq, file="hla_ref.RData")
