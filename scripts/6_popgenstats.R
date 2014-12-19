
load("5_refmapbias.rda")

# Directory to save plots
fig_dir <- "~/hla_tools/1kg/figs/"

library(reshape2)
library(ggplot2)

#-------------------------------------------------------------------------------
# 1. Calculate heterozygosity at 1000G (filtered SNPs vs. kept SNPs)

# Vector of allele frequencies
Pal_kg_l <- list()
Pal_pag_l <- list()
filt_l <- list()

for (l in loci){
    Pal_kg_l[[l]] <- c(Pal_kg_l[[l]],
                     reffreq[[l]]["global", polim.kg[[l]], "1000G"])
    Pal_pag_l[[l]] <- c(Pal_pag_l[[l]],
                     reffreq[[l]]["global", polim.kg[[l]], "PAG2014"])
    filt_l[[l]] <- c(filt_l[[l]], polim.kg[[l]] %in% badfreq[[l]])
}

Pal_kg <- unlist(Pal_kg_l)
Pal_pag <- unlist(Pal_pag_l)
filt <- unlist(filt_l)
filt2 <- filt
filt2[filt] <- "excluded"
filt2[!filt] <- "kept"

# Dataframe to be used for plotting
stat <- data.frame(Dat = c(rep("1000G", length(Pal_kg)),
                           rep("PAG2014", length(Pal_pag))),
                   Locus = c(rep(names(Pal_kg_l), sapply(Pal_kg_l, length)),
                             rep(names(Pal_pag_l), sapply(Pal_pag_l, length))),
                   Freq = c(Pal_kg, Pal_pag),
                   H = c(2 * Pal_kg * (1 - Pal_kg),
                         2 * Pal_pag * (1 - Pal_pag)),
                   Filt = rep(filt2,2))

#-------------------------------------------------------------------------------
# Plot 1000G Heterozygosity distribution for filtered and unfiltered sites

cbPalette <- c("#000000", "#D55E00", "#0072B2", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#CC79A7")

pdf(paste0(fig_dir,"H_denshist_kg.pdf"))
ggplot(subset(stat, Dat == "1000G"), aes(x = H)) +
geom_histogram(aes(y = ..density.., fill = Filt, colour = Filt),
               position = position_dodge(width = 0.025), binwidth = 0.05,
               alpha= 0.4) +
geom_histogram(aes(y = ..density.., fill = "all", colour = "all"),
               binwidth = 0.05, alpha = 0.2) +
scale_fill_manual(values=c("all" = cbPalette[1], "excluded" = cbPalette[2],
                           "kept" = cbPalette[3] )) +
scale_colour_manual(values=c("all" = NA, "excluded" = cbPalette[2],
                             "kept" = cbPalette[3] )) +
guides(colour = guide_legend(override.aes = list(alpha = c(0.2, 0.4, 0.4)))) +
ggtitle("Heterozygosity in 1000G") +
theme_classic() +
theme(legend.title=element_blank())
dev.off()

#-------------------------------------------------------------------------------
# Plot PAG2014 Heterozygosity distribution for filtered and unfiltered sites

pdf(paste0(fig_dir,"H_denshist_pag.pdf"))
ggplot(subset(stat, Dat == "PAG2014"), aes(x = H)) +
geom_histogram(aes(y = ..density.., fill = Filt, colour = Filt),
               position = position_dodge(width = 0.025), binwidth = 0.05, alpha = 0.4) +
geom_histogram(aes(y = ..density.., fill = "all", colour = "all"), binwidth = 0.05, alpha = 0.2) +
scale_fill_manual(values=c("all" = cbPalette[1], "excluded" = cbPalette[2],
                           "kept" = cbPalette[3] )) +
scale_colour_manual(values=c("all" = NA, "excluded" = cbPalette[2], "kept" = cbPalette[3] )) +
guides(colour = guide_legend(override.aes = list(alpha = c(0.2, 0.4, 0.4)))) +
ggtitle("Heterozygosity in PAG2014") +
theme_classic() +
theme(legend.title=element_blank())
dev.off()

#-------------------------------------------------------------------------------
# 2. Plot Heterozygosity distribution for filtered sites comparing 1000G vs. PAG

pdf(paste0(fig_dir,"H_denshist_filtered.pdf"))
ggplot(subset(stat, Filt == "excluded") , aes(x = H)) +
geom_histogram(aes(y = ..density.., fill = Dat, colour = Dat),
               position = position_dodge(width = 0.025), binwidth = 0.05,
               alpha= 0.4)+
ggtitle("Heterozygosity in excluded sites") +
theme_bw() +
theme(legend.title=element_blank())
dev.off()

#-------------------------------------------------------------------------------
# Calculate nucleotide diversity (pi) per locus

# number of sequences
n <- 930*2

# calculate pi per locus
# f is a vector of allele frequencies
# l is the sequence length
# n is the number of sequences
calc_nucdiv <- function(f, l, n){
    sum(2*f*(1-f))*(n/(n-1)) / l
}

nucdiv <- mapply(Pal_pag_l, sapply(reffreq, ncol), FUN = calc_nucdiv, n)

# Dataframe to be used for plotting
freqdevdf <- data.frame(
                        Deviation = unlist(sapply(freqdev, function(x)
                                                  x["global", ])),
                        NucDiv = rep(nucdiv, sapply(freqdev, ncol)))

pdf(paste0(fig_dir, "freqdev_boxplot.pdf"))
ggplot(freqdevdf, aes(NucDiv, abs(Deviation))) +
geom_boxplot(aes(fill = Locus))
dev.off()

stat$FreqDev <- unlist(sapply(freqdev, function(x) x["global", ]))

pdf(paste0(fig_dir, "freqdev_h_persite.pdf"), height = 6, width = 8)
ggplot(subset(stat, Dat == "PAG2014"), aes(x = H, y = FreqDev))+
geom_point() +
geom_smooth(method = lm, se = FALSE, colour = cbPalette[1]) +
labs(x = "H", y = "FE") +
theme_classic()
dev.off()

pdf(paste0(fig_dir, "freqdevabs_h_persite.pdf"), height = 6, width = 8)
ggplot(subset(stat, Dat == "PAG2014"), aes(x = H, y = abs(FreqDev)))+
geom_point() +
geom_smooth(method = lm, se = FALSE, colour = cbPalette[1]) +
labs(x = "H", y = "|FE|") +
theme_classic()
dev.off()

cor.test(stat[stat$Dat == "PAG2014", "H"],
         stat[stat$Dat == "PAG2014", "FreqDev"])
