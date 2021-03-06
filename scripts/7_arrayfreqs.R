###################################################################################
#
# Comparing HLA allele frequencies between Sanger sequencing and array genotyping
#
# Sep/2014
#
###################################################################################

source("3_af_accuracy.R")

#----------------------------------------------------------
# 1. Pierre data

data_dir <- "~/hla_tools/1kg/data/"
pierre.orig <- read.table(paste0(data_dir, "mhc.tab"), sep = "\t",
		   header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# Clean sample names
pierre.orig$Subject <- gsub("_.", "", pierre.orig$Subject)

# correct allele name 14XX
pierre.orig[pierre.orig=="14XX"] <- "14"

# Get list of individuals in Pierre data to filter Axiom data
#write.table(pierre$id,'../data/pierre_ind.txt',quote=F,row.names=F,col.names=F)

# Then, filter Pierre data for individuals in Axiom data
ind_ovlp <- read.table('~/Dropbox/Private/Mestrado/Analises/validate_hla1000g/data/ind_overlap_axiom_pierre.txt',
                       as.is=T,)
pierre.orig <- pierre.orig[pierre.orig$Subject %in%
                           ind_ovlp[,1],]

# Load consensus sequences for HLA alleles at Pierre's resolution generated by Vitor's script pierre.R
load("cons_seq_pierre.RDa")

# Recover HLA sequences from Pierre allele names
pierre <- list()
pierre$A <- sapply( pierre.orig[,5:6], function(x) cons_seq_pierre_A[x,],
                   simplify="array")
rownames(pierre$A) <- pierre.orig$Subject
pierre$B <- sapply( pierre.orig[,7:8], function(x) cons_seq_pierre_B[x,],
                   simplify="array")
rownames(pierre$B) <- pierre.orig$Subject
pierre$C <- sapply( pierre.orig[,9:10], function(x) cons_seq_pierre_C[x,] , simplify="array")
rownames(pierre$C) <- pierre.orig$Subject
pierre$DQB1 <- sapply( pierre.orig[,13:14], function(x) cons_seq_pierre_DQB[x,] , simplify="array")
rownames(pierre$DQB1) <- pierre.orig$Subject
pierre$DRB1 <- sapply( pierre.orig[,11:12], function(x) cons_seq_pierre_DRB[x,] , simplify="array")
rownames(pierre$DRB1) <- pierre.orig$Subject

# HLA-B -C -DQB1 and -DRB1 are on the minus strand. To allow comparison with
# the Axiom array data, convert nucleotides in the IMGT sequences, which by convention
# report the coding sequence, to the + strand

pierre[c("B", "C", "DQB1", "DRB1")] <- lapply(pierre[c("B", "C", "DQB1", "DRB1")],
    function(x){ chartr(new="ATCG", old="TAGC", x)-> x }
)

#-----------------------------------------------------------
# 2. Axiom data

axiom <- read.table('~/Dropbox/Private/Mestrado/Analises/validate_hla1000g//data/1000gSNPs.axiom.extMHC.pierre_ind.freq.txt',header=T,as.is=T)

axiom$GENE <- gsub("_ex.","",axiom$HLA_exon)
axiom$GENE <- gsub("HLA.","",axiom$GENE)
axiom$HLA_exon <- gsub(".+_ex","",axiom$HLA_exon)
axiom$REF <- sapply(strsplit(axiom$ALLELE,":"),"[",1)
axiom$REF_FREQ <- as.numeric(sapply(strsplit(axiom$ALLELE,":"),"[",2))
axiom$ALT<- sapply(strsplit(axiom$FREQ,":"),"[",1)
axiom$ALT_FREQ <- as.numeric(sapply(strsplit(axiom$FREQ,":"),"[",2))

#-----------------------------------------------------------
# 3. Select positions from axiom in Pierre data

# Coordinates of ARS exons in hg19. Note that some HLA loci are on the - strand
exars_hg19 <- list()
exars_hg19$A <- c(29910534:29910803,29911045:29911320)
exars_hg19$B <- c(31324734:31324465,31324219:31323944)
exars_hg19$C <- c(31239645:31239376,31239125:31238850)
exars_hg19$DQB1 <- c(32632844:32632575)
exars_hg19$DRB1 <- c(32552155:32551886)

# Relative positions of SNPs in Axiom data per locus
# (Relative to exons 2&3 coordinates)
loci <- c("A","B","C","DQB1","DRB1")
pos <- list()
for(l in loci){
    pos[[l]] <- as.numeric(na.omit(match(axiom$POS, exars_hg19[[l]])))
}

# Filter Pierre data keeping only positions present in Axiom
pier.new <- mapply(pierre,pos,FUN=function(x,y){x[,y,]})

#-----------------------------------------------------------
# 4. Allele frequencies in Pierre data

# List of individuals in each population
pop <- list()
for (p in unique(pierre.orig$Population)){
    pshort <- substr(p, 16, 18)
    pop[[pshort]] <- pierre.orig$Subject[pierre.orig$Population == p]
}
pop <- c(list(global = pierre.orig$Subject), pop)

#Load RData with reference sequences for ARS exons
load("hla_ref.RData")

# Organize ref sequences with only the polymorphic positions in Axiom data and 
# on + strand
rs <- list()
for(l in loci){
    rs[[l]] <- refseq[[l]][,pos[[l]]]
    if(l != "A"){
        rs[[l]] <- chartr("ATCG", "TAGC", rs[[l]])
    }
}

#Function to calculate reference allele frequency at each site for a given population
calcreffreq <- function(sangl,refl,pop){
    # rbinds the 2 alleles at each locus generating an array of
    #(2*nind_pop) lines by (axiom_sites) columns
    sang <- rbind(sangl[pop,,1],sangl[pop,,2])
    # vector to store frequencies at each site
    fr <- numeric(length=ncol(sang))
    names(fr) <- list(colnames(sang))
    # calculate frequency of reference allele at each site
    for(i in 1:ncol(sang)){
        fr[i] <- sum(sang[,i]==refl[i])/nrow(sang)
    }
    return(fr)
}

# Calculate reference allele frequency at each site per population
reffreq_axiom <- list()
for (l in loci){
    reffreq_axiom[[l]] <- matrix(nrow=length(pop),ncol=ncol(pier.new[[l]])
                           ,dimnames=list(names(pop),colnames(pier.new[[l]])))
    for(p in names(pop)){
        reffreq_axiom[[l]][p,] <- t(calcreffreq(pier.new[[l]],rs[[l]],pop[[p]]))
    }
}

#-----------------------------------------------------------
# 5. Compare allele frequencies

freqdev_axiom <- list()
for (l in loci){
    freqdev_axiom [[l]] <- axiom[axiom$GENE==l,"REF_FREQ"]-reffreq_axiom[[l]]["global",]
}

################################################################################

save.image("7_arrayfreqs.rda")

################################################################################
#-----------------------------------------------------------
# 6. Compare frequency deviation in Axiom and 1000G data (relative to Sanger)

# load rdata from kg validation
load("1_data_preproc.rda")

# select only sites present in both Axiom and 1000G
site_ovlp <- mapply(freqdev_axiom,freqdev,FUN=function(x,y) intersect(names(x),colnames(y)))

tables_frdev <- list()
for(l in loci){
    tables_frdev[[l]] <- round(rbind(kg=freqdev[[l]][1,site_ovlp[[l]]]
                                     ,axiom=freqdev_axiom[[l]][site_ovlp[[l]]])
                                ,6)
}

#-----------------------------------------------------------
# 7. Plot allele frequencies

# Directory to save plots
fig_dir <- "~/hla_tools/1kg/figs/"

freqthr <- 0.1
p='global'
ssres <- numeric()
rmse <- numeric()
for(l in loci){
    #Calculate global root-mean-square-error
    ssres[l] <- sum((reffreq_axiom[[l]][1,]-axiom[axiom$GENE==l,"REF_FREQ"])^2)
    rmse[l] <- sqrt(ssres[l]/ncol(reffreq_axiom[[l]]))
    #Plot globl allele freq per si per locus
    pdf(paste(fig_dir,"axiom_ref",l,p,"abs.pdf",sep="_"),useDingbats=F)
    par(mar=c(5,5,5,7),xpd=T)
    plot(x=reffreq_axiom[[l]][p,]
         ,y=axiom[axiom$GENE==l,"REF_FREQ"],
         ylab="REF freq Axiom", xlab="REF freq HLAdat",ylim=c(0,1.1),xlim=c(0,1.1),
         main=paste0("HLA-",l," ",p),
         xaxs="i",yaxs='i',bty='n')
    lines(0:1,0:1,col="dimgray")
    text(y=axiom[axiom$GENE==l,"REF_FREQ"]+0.02,
         x=reffreq_axiom[[l]][p,],
         labels=colnames(reffreq_axiom[[l]]), cex=0.5)
    if(length(site_ovlp[[l]]>0)){
        points(x=reffreq[[l]][p,site_ovlp[[l]],'sanger']
             ,y=reffreq[[l]][p,site_ovlp[[l]],'kg'],col='gray50')
        text(y=reffreq[[l]][p,site_ovlp[[l]],'kg']+0.02,
             x=reffreq[[l]][p,site_ovlp[[l]],'sanger'],
             labels=site_ovlp[[l]], cex=0.5)
    }
    mtext(paste0("RMSE =",round(rmse[l],digits=2)),cex=.7)
    lines(c(0,1-freqthr),c(0+freqthr,1)
           ,type='l',lty=2
           ,col='steelblue'
           )
    lines(c(0+freqthr,1),c(0,1-freqthr)
           ,type='l',lty=2
           ,col='steelblue'
           )
    dev.off()
}

# sites present in both Axiom and 1000G
site_ovlp1 <- mapply(freqdev,freqdev_axiom,FUN=function(x,y) colnames(x) %in% names(y))
site_ovlp2 <- mapply(freqdev,freqdev_axiom,FUN=function(x,y) names(y) %in% colnames(x))

dat <- rbind(data.frame(Site = unlist(sapply(freqdev_axiom, names)),
                        Site_ovl = unlist(site_ovlp2),
                        FreqDev=unlist(freqdev_axiom), Source="Axiom"),
             cbind(Site = unlist(sapply(freqdev, colnames)),
                   Site_ovl = unlist(site_ovlp1),
                   FreqDev=unlist(sapply(freqdev, function(x) x[1,])), Source = "1000 Genomes"))


library("ggplot2")
pdf(paste0(fig_dir, "axiom_freqdev1.pdf"))
ggplot(dat, aes(factor(Source), as.numeric(FreqDev))) +
geom_jitter(aes(colour = Site_ovl), position = position_jitter(width = .3), alpha = 0.6) +
scale_colour_manual(values = c("black", "red")) +
theme_classic()
dev.off()

pdf(paste0(fig_dir, "axiom_freqdev2.pdf"))
ggplot(dat[-15,], aes(factor(Source), as.numeric(FreqDev))) +
geom_hline(yintercept = 0, colour = "gray") +
geom_jitter(aes(colour = Site_ovl), position = position_jitter(width = .3), alpha = 0.6) +
scale_colour_manual(values = c("black", "red"),
                    breaks = c(FALSE, TRUE),
                    labels = c("Axiom OR 1000G", "Axiom AND 1000G"),
                    name = "SNP source") +
labs(x = "Dataset", y = "Frequency difference to PAG2014 (FE)") +
theme_classic()
dev.off()

test_freqdif <- permtest(abs(as.numeric(dat$FreqDev[-15])),
dat$Source[-15] == "Axiom")
