library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(data.table)
library(ggiraphExtra)
library(webr)
library(plotly)
require(ggiraph)
require(plyr)
library(orca)


setwd("20260528_NUR77_paper_codes")

writeLines(capture.output(sessionInfo()), "sessionInfo_chipseeker_annotation.txt")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

###########################################
#For overlapping C3,C5,C6 CsNB samples
###########################################
dcsnb <- data.frame(fread(file="20260528_NUR77_paper_codes/C3C5C6_CsNB_overlapped_peaks_final.bed"))

#############################################################################################
#Extract peak sites that have 1,2, or 3 samples
###########################################################################################
dcsnb_all <- dcsnb[, c("chrom","start","end")]
dcsnb_all$seq <- paste0("peak_",seq(1:nrow(dcsnb_all)))
fwrite(dcsnb_all, file="C3C5C6_CsNB_overlapped_trimmed.bed",sep = "\t", col.names = FALSE)
peak_all <- readPeakFile("20260528_NUR77_paper_codes/C3C5C6_CsNB_overlapped_trimmed.bed")

peakAnno_all <- annotatePeak(peak_all, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")

#Visualize genomic annotation
annoDetails_all <- data.frame(peakAnno_all@annoStat)
annoDetails_all$category <- c("promoter","promoter","promoter","body","body","body","body","body","body","intergenic","intergenic")


#pie chart
plot_all <- plot_ly(annoDetails_all) %>%
  add_pie(annoDetails_all, labels = ~`Feature`, values = ~`Frequency`,
          domain = list(
            x = c(0.15, 0.85),
            y = c(0.15, 0.85)),
          sort = F) %>% layout(title = 'C3C5C6_overlapped_peaks(one or more of C3/C5/C6)')
plot_all

#Extract peak sites that have atleast 2 samples
dcsnb_2 <- dcsnb[dcsnb$num==2, c("chrom","start","end")]
dcsnb_2$seq <- paste0("peak_",seq(1:nrow(dcsnb_2)))
fwrite(dcsnb_2, file="C3C5C6_CsNB_overlapped2_trimmed.bed",sep = "\t", col.names = FALSE)
peak_2 <- readPeakFile("20260528_NUR77_paper_codes/C3C5C6_CsNB_overlapped2_trimmed.bed")

peakAnno_2 <- annotatePeak(peak_2, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")

#Visualize genomic annotation
annoDetails_2 <- data.frame(peakAnno_2@annoStat)
annoDetails_2$category <- c("promoter","promoter","promoter","body","body","body","body","body","body","intergenic","intergenic")


#pie chart
plot_2 <- plot_ly(annoDetails_2) %>%
  add_pie(annoDetails_2, labels = ~`Feature`, values = ~`Frequency`,
          domain = list(
            x = c(0.15, 0.85),
            y = c(0.15, 0.85)),
          sort = F) %>% layout(title = 'C3C5C6_overlapped_peaks(2 of C3/C5/C6)')
plot_2


#Extract peak sites that have all 3 samples
dcsnb_3 <- dcsnb[dcsnb$num==3, c("chrom","start","end")]


dcsnb_3$seq <- paste0("peak_",seq(1:nrow(dcsnb_3)))
fwrite(dcsnb_3, file="C3C5C6_CsNB_overlapped3_trimmed.bed",sep = "\t", col.names = FALSE)
peak_3 <- readPeakFile("20260528_NUR77_paper_codes/C3C5C6_CsNB_overlapped3_trimmed.bed")

peakAnno_3 <- annotatePeak(peak_3, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Hs.eg.db")

#Visualize genomic annotation
annoDetails_3 <- data.frame(peakAnno_3@annoStat)
annoDetails_3$category <- c("promoter","promoter","promoter","body","body","body","body","body","body","intergenic","intergenic")


#pie chart
plot_3 <- plot_ly(annoDetails_3) %>%
  add_pie(annoDetails_3, labels = ~`Feature`, values = ~`Frequency`,
          domain = list(
            x = c(0.2, 0.8),
            y = c(0.2, 0.8)),
          sort = F) %>% layout(title = 'C3C5C6_overlapped_peaks(3 of C3/C5/C6)')
plot_3


#Extract peak sites that have all 2 AND 3 samples
dcsnb_2_3 <- dcsnb[dcsnb$num==2 | dcsnb$num==3, c("chrom","start","end")]

dcsnb_2_3$seq <- paste0("peak_",seq(1:nrow(dcsnb_2_3)))
fwrite(dcsnb_2_3, file="C3C5C6_CsNB_overlapped2_3_trimmed.bed",sep = "\t", col.names = FALSE)
peak_2_3 <- readPeakFile("20260528_NUR77_paper_codes/C3C5C6_CsNB_overlapped2_3_trimmed.bed")

peakAnno_2_3 <- annotatePeak(peak_2_3, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Hs.eg.db")

#Visualize genomic annotation
annoDetails_2_3 <- data.frame(peakAnno_2_3@annoStat)
annoDetails_2_3$category <- c("promoter","promoter","promoter","body","body","body","body","body","body","intergenic","intergenic")


#pie chart
plot_2_3 <- plot_ly(annoDetails_2_3) %>%
  add_pie(annoDetails_2_3, labels = ~`Feature`, values = ~`Frequency`,
          domain = list(
            x = c(0.2, 0.8),
            y = c(0.2, 0.8)),
          sort = F) %>% layout(title = 'C3C5C6_overlapped_peaks(2 or 3 of C3/C5/C6)')
plot_2_3


###############################################
#For C3/C5/C6 UT samples
###############################################

dut<- data.frame(fread(file="20260528_NUR77_paper_codes/C3C5C6_UT_overlapped_peaks_final.bed"))

#Having all 3 samples overlapping
dut_3 <- dut[dut$num==3, c("chrom","start","end")]
dut_3$seq <- paste0("peak_",seq(1:nrow(dut_3)))
fwrite(dut_3, file="C3C5C6_UT_overlapped3_trimmed.bed",sep = "\t", col.names = FALSE)

#Having 2 and 3 samples overlapping
dut_2_3 <- dut[dut$num==2 | dut$num==3, c("chrom","start","end")]
dut_2_3$seq <- paste0("peak_",seq(1:nrow(dut_2_3)))
fwrite(dut_2_3, file="C3C5C6_UT_overlapped_2_3_trimmed.bed",sep = "\t", col.names = FALSE)


#############################################
#For differentially bound sites
#Do the same analysis as above on the differentially bound sites (dbObj_count_norm.DB.batch2.bed) and 
#only for the significant peaks from the differentially bound sites (dbObj_count_norm.DB.Significant.batch2.bed)
#############################################

diffpeaks <- readPeakFile("20260528_NUR77_paper_codes/dbObj_count_norm.DB.batch2.bed")

diffpeakAnno <- annotatePeak(diffpeaks, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")

#Visualize genomic annotation
annoDetails_diffpeaks <- data.frame(diffpeakAnno@annoStat)
annoDetails_diffpeaks$category <- c("promoter","promoter","promoter","body","body","body","body","body","body","intergenic","intergenic")


annoDetails_diffpeaks$legend_label <- paste(annoDetails_diffpeaks$Feature, paste0(round(annoDetails_diffpeaks$Frequency,2),"%"), sep=" ")

#pie chart
plot_diff <- plot_ly(annoDetails_diffpeaks) %>%
  add_pie(annoDetails_diffpeaks, labels = ~`legend_label`, values = ~`Frequency`,
          domain = list(
            x = c(0.15, 0.85),
            y = c(0.15, 0.85)),
          sort = F,textinfo = "none") %>% layout(title = 'Differntially bound peaks')
plot_diff


#Significantly differentially bound sites
Sigdiffpeaks <- readPeakFile("20260528_NUR77_paper_codes/dbObj_count_norm.DB.Significant.batch2.bed")

SigdiffpeakAnno <- annotatePeak(Sigdiffpeaks, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")

#Visualize genomic annotation
annoDetails_sigdiffpeaks <- data.frame(SigdiffpeakAnno@annoStat)
annoDetails_sigdiffpeaks$category <- c("promoter","promoter","promoter","body","body","body","body","body","body","intergenic")

#pie chart
plot_sigdiff <- plot_ly(annoDetails_sigdiffpeaks) %>%
  add_pie(annoDetails_sigdiffpeaks, labels = ~`Feature`, values = ~`Frequency`,
          domain = list(
            x = c(0.2, 0.8),
            y = c(0.2, 0.8)),
          sort = F) %>% layout(title = 'Significantly differentially bound peaks')
plot_sigdiff

annoDetails_sigdiffpeaks$legend_label <- paste(annoDetails_sigdiffpeaks$Feature, paste0(round(annoDetails_sigdiffpeaks$Frequency,2),"%"), sep=" ")
plot_sigdiff <- plot_ly(annoDetails_sigdiffpeaks) %>%
  add_pie(annoDetails_sigdiffpeaks, labels = ~`legend_label`, values = ~`Frequency`,
          domain = list(
            x = c(0.15, 0.85),
            y = c(0.15, 0.85)),
          sort = F, textinfo = "none") %>% layout(title = 'Significantly differentially bound peaks')

plot_sigdiff <- plot_ly(annoDetails_sigdiffpeaks) %>%
  add_pie(annoDetails_sigdiffpeaks, labels = ~`legend_label`, values = ~`Frequency`,
          domain = list(
            x = c(0.15, 0.85),
            y = c(0.15, 0.85)),
          sort = F) %>% layout(title = 'Significantly differentially bound peaks')
plot_sigdiff

##################################################
#Plot heatmap profiles of peaks
#################################################
#Peak count for all the samples

# Define a window extending 3kb upstream and 3kb downstream of the TSS
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

tagMatrixList <- list()
files <- list()
files[[1]] <- readPeakFile("/data_c2/data_n3/home/menaka/Projects/20240612_ShiYong_chipseq_batch2/paper_submission_related_coded_data_visuals/C3_CsNB_out_peaks.narrowPeak.bed6.bed")
peak1 <- getTagMatrix(files[[1]], windows=promoter)
files[[2]] <- readPeakFile("/data_c2/data_n3/home/menaka/Projects/20240612_ShiYong_chipseq_batch2/paper_submission_related_coded_data_visuals/C5_CsNB_out_peaks.narrowPeak.bed6.bed")
peak2 <- getTagMatrix(files[[2]], windows=promoter)
files[[3]] <- readPeakFile("/data_c2/data_n3/home/menaka/Projects/20240612_ShiYong_chipseq_batch2/paper_submission_related_coded_data_visuals/C6_CsNB_out_peaks.narrowPeak.bed6.bed")
peak3 <- getTagMatrix(files[[3]], windows=promoter)

files[[4]]<- readPeakFile("/data_c2/data_n3/home/menaka/Projects/20240612_ShiYong_chipseq_batch2/paper_submission_related_coded_data_visuals/C3_UT_out_peaks.narrowPeak.bed6.bed")
peak4 <- getTagMatrix(files[[4]], windows=promoter)
files[[5]] <- readPeakFile("/data_c2/data_n3/home/menaka/Projects/20240612_ShiYong_chipseq_batch2/paper_submission_related_coded_data_visuals/C5_UT_out_peaks.narrowPeak.bed6.bed")
peak5 <- getTagMatrix(files[[5]], windows=promoter)
files[[6]] <- readPeakFile("/data_c2/data_n3/home/menaka/Projects/20240612_ShiYong_chipseq_batch2/paper_submission_related_coded_data_visuals/C6_UT_out_peaks.narrowPeak.bed6.bed")
peak6 <- getTagMatrix(files[[6]], windows=promoter)

tagMatrixList <- list("C3_CsNB"=peak1, "C5_CsNB"=peak2, "C6_CsNB"=peak3,
                      "C3_UT"=peak4,"C5_UT"=peak5,"C6_UT"=peak6)


#Peak heatmaps - make sure you have the ChIPseeker_1.40.0 version for the below to produce the profile Heatmap!
pdf("profiles_of_peaks_heatmaps.pdf", width=18, height=10)
tagHeatmap(tagMatrixList, palette = "YlGn") #"YlGn",RdBu"
dev.off()



