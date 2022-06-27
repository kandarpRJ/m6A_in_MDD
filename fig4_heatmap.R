library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("~/data/MDD_review/metaanalysis")

diffpeaks<-read.csv("temporal/DiffMod.csv")
peakcounts<-read.csv("temporal/ADDInfo/ADDInfo_ReadsCount.csv")
selDEres<-read.table("Suppl2.xls", sep = "\t", header = T)

peaks_df<-merge(peakcounts, diffpeaks, by.x="X", by.y="name")
peaks_df<-peaks_df[,c("X", "geneID", "human_young_FrontalCortex_1.m6A.bam", "human_young_FrontalCortex_2.m6A.bam", 
                      "human_young_FrontalCortex_1.input.bam", "human_young_FrontalCortex_2.input.bam", 
                      "human_old_FrontalCortex_1.m6A.bam", "human_old_FrontalCortex_2.m6A.bam",
                      "human_old_FrontalCortex_1.input.bam", "human_old_FrontalCortex_2.input.bam",
                      "DiffModLog2FC", "pvalue", "padj")]
peaks_df$geneID<-gsub("\\..*", "", peaks_df$geneID)

gls<-bitr(peaks_df$geneID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db, drop = FALSE)

peaks_df<-merge(peaks_df, gls, by.x="geneID", by.y="ENSEMBL", all.x=T)

keep<-peaks_df$SYMBOL %in% rownames(selDEres[1:100,])

dis_peaks<-peaks_df[keep,]
colnames(dis_peaks)<-c("geneID", "peak", "YoungAdult1_IP", "YoungAdult2_IP", "YoungAdult1_INPUT", "YoungAdult2_INPUT", 
                       "OldAdult1_IP", "OldAdult2_IP", "OldAdult1_INPUT", "OldAdult2_INPUT", "Log2FC", "pvalue", "padj", "SYMBOL")
forcol<-dis_peaks$pvalue
forcol[forcol>0.05]<-1
#padj_col = colorRamp2(sort(dis_peaks$pvalue, decreasing = T), hcl.colors(length(dis_peaks$pvalue),"PuBu"))
pcol<-colorRamp2(c(0,0.05,0.08,0.1), c("purple", "white", "orange", "darkorange"))
row_ha<-rowAnnotation(Pvalue=dis_peaks$pvalue, col=list(Pvalue=pcol))

heatmapcol<-colorRamp2(c(0,5,10), c("darkred", "white", "darkblue"))

tiff("DiffMod_IPCounts_heatmap_top100DE.tif", res = 300, width = 3000, height = 3200)
Heatmap(log(dis_peaks[,c(3,4,7,8)]), name="log(IP_counts)", right_annotation = row_ha, show_row_names = F, 
        cluster_rows = T, cluster_columns = F, row_names_gp = gpar(fontsize=14), column_names_gp = gpar(fontsize=14))
dev.off()



diffgenes<-read.table("old_vs_young.input.de.res.txt", sep = "\t", row.names = NULL, header = T)
genecounts<-read.csv("temporal/genes.csv", row.names = NULL)
genecounts$gene_id<-gsub("ENSG\\d*\\.\\d*\\|","",genecounts$gene_id)
gcmat<-merge(diffgenes, genecounts, by.x="rowname", by.y="gene_id")

selgc<-gcmat[gcmat$rowname %in% selDEres$Row.names[1:100],] 
rownames(selgc)<-selgc$rowname
selgc<-selgc[,-1]
selgc<-selgc[selDEres$Row.names[1:100],]
selgc$significance<-"No"
selgc$significance[selgc$P.Value<0.05]<-"Yes"

rha<-rowAnnotation(Pvalue=selgc$P.Value, col=list(Pvalue=pcol))

tiff("DiffMod_InputCounts_heatmap_top100DE.tif", res = 300, width = 2500, height = 4000)
Heatmap(log(selgc[,c("young_1", "young_2", "old_1", "old_2")]), name="log(Input_counts)", right_annotation = rha, 
        cluster_columns = F, cluster_rows = F, row_names_gp = gpar(fontsize=14), column_names_gp = gpar(fontsize=14),
        column_labels = c("YoungAdult1_Input", "YoungAdult2_Input", "OldAdult1_Input", "OldAdult2_Input"))
dev.off()
