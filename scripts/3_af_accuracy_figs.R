################################################################################
#
# Figs for 3_af_accuracy.R
#
# Allele frequency deviation & reference allele bias
#
################################################################################
#-------------------------------------------------------------------------------
# 0. Load data
source("3_af_accuracy.R")

# Directory to save plots
fig_dir <- "../figures/"

# Colorblind safe palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
               "#D55E00", "#CC79A7")

# Wordcloud package to avoid overplotting
install.packages("wordcloud")
library(wordcloud)
source("mytextplot.R")

#-------------------------------------------------------------------------------
# 1. Plot scatterplot of frequency at 1000G vs. PAG2014

freqthr <- 0.1

pop_codes <-c("Global", "ASW", "GBR","CHB+JPT","CLM", "FIN", "CHS", "TSI", 
              "LWK", "MXL", "CEU", "PUR", "YRI")

for(l in loci){
    #Plot allele freq per site per locus per pop
    for(p in 1:length(pop)){
        p_name <- names(pop)[p]
        p_code <- pop_codes[p]
        pdf(paste0(fig_dir, "ref", l, p_name, ".pdf"), useDingbats = F)
        par(mar = c(4, 4, 4, 0), xpd = T)
        plot(x = reffreq[[l]][p_name, , "PAG2014"],
             y = reffreq[[l]][p_name, , "1000G"],
             xlab = "REF freq PAG2014", ylab = "REF freq 1000G",
             xlim = c(0, 1.1), ylim = c(0, 1.1), xaxs = "i", yaxs = 'i',
             main = paste0("HLA-", l, " ", p_code), bty = 'n', pch = 16,
             col = "gray65")
        lines(0:1, 0:1, col = "gray")
        mytextplot(x = reffreq[[l]][p_name, polim.kg[[l]], "PAG2014"],
             y = reffreq[[l]][p_name, polim.kg[[l]], "1000G"],
             words = polim.kg[[l]], cex = 0.8, new = F)
        mtext(paste0("MAE=", round(mae[p_name,l], digits = 2)))
        points(c(0, 1 - freqthr), c(0 + freqthr, 1), type = 'l', lty = 2,
               col = cbPalette[6])
        points(c(0 + freqthr, 1), c(0, 1 - freqthr), type = 'l', lty = 2,
               col = cbPalette[6])
        dev.off()
    }
}

#-------------------------------------------------------------------------------
# 2. Bland-Altman plot 1000G vs. PAG2014

for(l in loci){
    #Plot allele freq per site per locus per pop
    for(p in 1:length(pop)){
        p_name <- names(pop)[p]
        p_code <- pop_codes[p]
        pdf(paste0(fig_dir, "BlAlt_", l, p_name, ".pdf"), useDingbats = F)
        par(mar = c(4, 4, 4, 0), xpd = F)
        f1 <- reffreq[[l]][p_name,polim.kg[[l]] , "1000G"]
        f2 <- reffreq[[l]][p_name,polim.kg[[l]] , "PAG2014"]
        diff <- f1-f2
        mdiff <- mean(diff)
        sddiff <- sd(diff)
        plot(x = (f1+f2)/2, y = diff,
             xlab = "Average frequency", ylab = "Frequency difference",
             xlim = c(0,1.15), ylim = c(-0.6, 0.6), xaxs = "i", yaxs = "i",
             main = paste0("HLA-", l, " ", p_code),
             bty = 'n', pch = 16, col = "gray65")
        lines(c(0,1), c(0,0), col = cbPalette[7])
        lines(c(0,1), c(mdiff, mdiff), col = "gray")
        lines(c(0,1), rep(mdiff + 1.96 * sddiff, 2), lty = 2, col = "gray")
        lines(c(0,1), rep(mdiff - 1.96 * sddiff, 2), lty = 2, col = "gray")
        mytextplot(x = (f1 + f2)/2, y = diff,
             words = polim.kg[[l]], cex = 0.8, new = F)
        dev.off()
    }
}

#----------------------------------------------------------------------------------
# 3. Comparing mismatch and frequency deviation

for (l in loci){
    pdf(paste0(fig_dir, "mism_freqdev_", l, ".pdf"), height = 6, width = 10)
    plot(apply(mism[[l]], 2, sum)/nrow(mism[[l]]), abs(freqdev[[l]]["global",]),
         main = paste0("HLA-", l), xlab = "Genotype mismatches per site",
         ylab = "Mean frequency difference (MAE)", 
         xlim = c(0,1), ylim = c(0,.6), pch = 16, col = "grey65", bty="l")
    mytextplot(x = apply(mism[[l]], 2, sum)/nrow(mism[[l]]),
               y = abs(freqdev[[l]]["global",]),
               words = colnames(mism[[l]]), new = F, cex = .8,
               xlim = c(0,1), ylim = c(0,.6))
    dev.off()
}

# correlation:
pmism <- numeric()
fd <- numeric()
for (l in loci){
    pmism <- c(pmism, apply(mism[[l]], 2, sum)/nrow(mism[[l]]))
    fd <- c(fd, abs(freqdev[[l]]["global",]))
}

cor.test(pmism, fd)
