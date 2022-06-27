library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

setwd("~/data/MDD_review/metaanalysis")

meta<-read.csv("MDD_study_metadata.csv")
selDEres<-read.table("Suppl2.xls", sep = "\t", header = T)

col_fun1 = colorRamp2(selDEres$adjP.fishercomb[1:100], hcl.colors(length(selDEres$adjP.fishercomb[1:100]),"PuRd"))
col_fun2 = colorRamp2(selDEres$adjP.invnormcomb[1:100], hcl.colors(length(selDEres$adjP.invnormcomb[1:100]),"PuBu"))
#col_fun3 = colorRamp2(selDEres$m6A_count[1:100], hcl.colors(length(selDEres$m6A_count[1:100]),"heat"))
col_fun3<-colorRamp2(c(0,1,20), c("purple", "yellow", "red"))

brcol<-brewer.pal(n=8, name = "Paired")
names(brcol)<-unique(meta$Brain.region[meta$file %in% paste(colnames(selDEres[,2:10]), ".txt", sep="")])

sbrcol<-brewer.pal(n=8, name = "Set3")
names(sbrcol)<-unique(meta$Subregion[meta$file %in% paste(colnames(selDEres[,2:10]), ".txt", sep="")])
sbrcol[4]<-"black"

row_ha<-rowAnnotation(Pval_FishComb=selDEres$adjP.fishercomb[1:100], Pval_InvNormComb=selDEres$adjP.invnormcomb[1:100], 
                      m6A_count=selDEres$m6A_count[1:100], col=list(fishcomb.Pval=col_fun1, invnormcomb.Pval=col_fun2, m6A_count=col_fun3),
                      annotation_legend_param=list(Pval_FishComb=list(labels_gp = gpar(fontsize = 14)),
                                                   Pval_InvNormComb=list(labels_gp = gpar(fontsize = 14)),
                                                   m6A_count=list(labels_gp = gpar(fontsize = 14))))

tiff("Combined_pval_heatmap.tif", width = 3000, height = 4700, res = 300)
h<-Heatmap(selDEres[1:100,2:10], name = "Significance (Yes-red/No-blue)", 
           #column_labels = gsub("_.*", "", colnames(selDEres[2:10])),
           column_labels = colnames(selDEres[2:10]),
           #top_annotation = ha, 
           right_annotation = row_ha, column_names_gp = grid::gpar(fontsize=14),
           row_names_gp = grid::gpar(fontsize=14), heatmap_legend_param = list(labels_gp = gpar(fontsize = 14)),
           show_row_names = T, cluster_rows = F, cluster_columns = F)
draw(h, annotation_legend_side="top")
dev.off()
