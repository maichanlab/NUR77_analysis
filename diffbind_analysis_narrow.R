library(DiffBind)
library(tidyverse)


setwd("20260528_NUR77_paper_codes")

writeLines(capture.output(sessionInfo()), "sessionInfo_diffbind_analysis.txt")

#Read in narrowpeaks pertaining to samples described in the samplesheet.csv
samples <- read.csv("20260528_NUR77_paper_codes/samplesheet.csv")


#Extract a consensus peakset with minOverlap=2
#Further analysis will be carried out on this consensus set 
dbObj <- dba(sampleSheet=samples) 

#Compute count information from alignment files, for each of 
#the peaks/regions in the consensus set. 
dbObj_count <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

dbObj_count
#6 Samples, 53313 sites 

#For each sample, multiplying the value in the Reads column by the corresponding FRiP value
#will yield the number of reads that overlap a consensus peak. 
info <- dba.show(dbObj_count)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes

#PCA plot
# pdf('pca_sample_clustering_for_consensusSites.pdf')
# dba.plotPCA(dbObj_count,  attributes=DBA_CONDITION)
# dev.off()


#correlation heatmap to evaluate the relationship between samples at consensus sites.
# pdf('heatmap_sample_correlation_consensusSites.pdf')
# plot(dbObj_count)
# dev.off()

dbObj_count_norm <- dba.normalize(dbObj_count)

###############################
#Differential Enrichment Analysis
###############################
#Contrast checking : treated vs untreated (treated as group1 and untreated as group2)
#Untreated as the baseline

dbObj_count_norm <- dba.contrast(dbObj_count_norm,reorderMeta=list(Condition="untreated"))

dbObj_count_norm <- dba.analyze(dbObj_count_norm, bBlacklist = FALSE, bGreylist = FALSE)

dba.show(dbObj_count_norm , bContrasts=TRUE)


#Plot heatmap
# pdf("heatmap_w_diffbound_sites.pdf")
# plot(dbObj_count_norm, contrast=1)
# dev.off()


#Retrieve the differentially bound sites, returned as a GRanges object
res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
dbObj_count_norm.DB <- dba.report(dbObj_count_norm)
dbObj_count_norm.DB

#2664 of the 53313 sites are identified as being significantly differentially bound (DB) sites using the default threshold of
#FDR <= 0.05 with fold=0 (log2(1)); ie  FC>1)


# Write to file
out <- as.data.frame(dbObj_count_norm.DB)
out$seq <- paste("seq",seq(1:nrow(out)), sep="-")
out$contrast <- "treated vs untreated"
write.table(out, file="treated_vs_untreated_deseq2.txt", sep="\t", quote=F, row.names=F)

out_sig <- out
out_sig$peak_status <- "Not Sig"
out_sig$peak_status[out_sig$Fold > 0.6 & out_sig$FDR < 0.05] <- "FDR Up"
out_sig$peak_status[out_sig$Fold < -0.6 & out_sig$FDR < 0.05] <- "FDR Down"
out_sig$peak_status[out_sig$Fold > 0.6 & out_sig$FDR < 0.05] <- "Sig Up"
out_sig$peak_status[out_sig$Fold < -0.6 & out_sig$FDR < 0.05] <- "Sig Down"
write.table(out_sig, file="treated_vs_untreated_deseq2_for_GREAT.txt", sep="\t", quote=F, row.names=F)


#The final two columns give
#confidence measures for identifying these sites as differentially bound, with a raw p-value and
#a multiple-testing corrected FDR in the final column (also calculated by the DESeq2 analysis)

#Create BED files for each set of significant regions identified by DESeq2, 
#separating them based on the gain or loss of enrichment. 
#Write these regions to file and use as input for downstream visualization.

#Create bed files for each keeping only significant peaks of Treated (FDR < 0.05 & +ve Fold (Fold>0))
treated_enrich <- out %>% 
  dplyr::filter(FDR < 0.05 & Fold > 0) %>% 
  dplyr::select(seqnames, start, end, seq)

# Write to file (2048)
write.table(treated_enrich, file="treated_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

#Create bed files for each keeping only significant peaks of Untreated (FDR < 0.05 & -ve Fold (Fold<0))
untreated_enrich <- out %>% 
  dplyr::filter(FDR < 0.05 & Fold < 0) %>% 
  dplyr::select(seqnames, start, end, seq)

# Write to file (616)
write.table(untreated_enrich, file="untreated_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

#Write the GRanges object to a bedfile, may need it for downstream processing.
#The only trick is remembering the BED uses 0-based coordinates.


dt <- data.frame(seqnames=seqnames(dbObj_count_norm.DB),
                 starts=start(dbObj_count_norm.DB)-1,
                 ends=end(dbObj_count_norm.DB))

dt$seq <- paste("seq",seq(1:nrow(dt)), sep="-")
write.table(dt, file="dbObj_count_norm.DB.bed", quote=F, sep="\t", row.names=F, col.names=F)


library(Repitools)
#Convert Grange object to a data.frame
dbObj_count_norm.DB.df <- annoGR2DF(dbObj_count_norm.DB)
dbObj_count_norm.DB.df$seq <- paste("seq",seq(1:nrow(dbObj_count_norm.DB.df)), sep="-")
write.csv(dbObj_count_norm.DB.df, file="dbObj_count_norm.DB.all.csv", quote=F)

#Compare the number of differential bound sites that have enriched ER binding in
#the untreated samples and those with enriched binding in the treated samples:

#ER binding in the treated samples
sum(dbObj_count_norm.DB$Fold>0) #2048

#ER binding in the untreated samples
sum(dbObj_count_norm.DB$Fold<0) #616

#PCA plots
#A PCA plot using only the 2664 differentially bound sites, using
#an FDR threshold of 0.05

# pdf("pca_plots_diff_bound_batch2.pdf")
# dba.plotPCA(dbObj_count_norm, attributes=DBA_CONDITION, contrast=1,  dotSize=2.5, bLog = TRUE)
# dev.off()


#Use labels to find outlier sites
sigSites <- dba.plotVolcano(dbObj_count_norm, fold=log2(1.5))
#Use these sigSites for GREAT
dsig <- data.frame(seqnames=seqnames(sigSites),
                 starts=start(sigSites)-1,
                 ends=end(sigSites))

dsig$seq <- paste("seq",seq(1:nrow(dsig)), sep="-")
write.table(dsig, file="dbObj_count_norm.DB.Significant.bed", quote=F, sep="\t", row.names=F, col.names=F)


#how are reads distributed amongst the different classes of differentially bound sites and
#sample groups
# pvals <- dba.plotBox(dbObj_count_norm)
# pvals


#Heatmap to view binding affinity directly in the differentially bound sites
library(viridis)
hmap <- viridis(25)
pdf("heatmap_diff_binding_sites_withscaling.pdf", width=6, height=10)
readscores <- dba.plotHeatmap(dbObj_count_norm, contrast=1, correlations=FALSE,scale="row", colScheme = hmap,cexCol = 1,
                              ColAttributes = c(DBA_CONDITION, DBA_REPLICATE), key=FALSE)
dev.off()


