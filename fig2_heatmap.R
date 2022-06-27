library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)

setwd("~/data/MDD_review/metaanalysis")

filelist = list.files(pattern = ".*.txt")

m6amachine<-read.table("m6a_machine")
study.meta<-read.csv("MDD_study_metadata.csv")

#assuming tab separated values with a header    
datalist = lapply(filelist, function(x) read.csv(x, sep = "\t", header=T, row.names = NULL)) 

datalist<-lapply(datalist, function (x) distinct(x, Gene.symbol, .keep_all = T))

#assuming the same header/columns for all files
datafr = Reduce(function (x, y) merge (x, y, by="Gene.symbol", all=T), datalist)

selfdr<-grep("adj.Pval|adj.P.Val", colnames(datafr))
selpval<-grep("P.Value", colnames(datafr))
selfc<-grep("logFC", colnames(datafr))

fdr.mat<-datafr[,c(1,selfdr)]
pval.mat<-datafr[,c(1, selpval)]
lfc.mat<-datafr[,c(1,selfc)]

filelist<-gsub(".txt", "", filelist)
filelist<-gsub("_.*", "", filelist)

colnames(fdr.mat)<-c("Gene.symbol", filelist)
colnames(pval.mat)<-c("Gene.symbol", filelist)
colnames(lfc.mat)<-c("Gene.symbol", filelist)

m6a.fdr<-fdr.mat[fdr.mat$Gene.symbol %in% m6amachine$V1,]
m6a.pval<-pval.mat[pval.mat$Gene.symbol %in% m6amachine$V1,]
m6a.lfc<-lfc.mat[lfc.mat$Gene.symbol %in% m6amachine$V1,]
plotmat<-m6a.lfc


for(i in 2:ncol(plotmat)) {
    plotmat[(m6a.pval[,i]>=0.05 | is.na(m6a.pval[,i]))==TRUE,i]<-NA
}

tmat<-t(plotmat[,2:33])
tmat<-apply(tmat,2,as.numeric)
colnames(tmat)<-plotmat$Gene.symbol
rownames(tmat)<-colnames(plotmat[,2:33])

tmat<-cbind(study.meta[order(study.meta$Study.series),], tmat)
tmat<-tmat[with (tmat, order(tmat$Brain.region, tmat$Subregion, tmat$Sex)),]

rseq<-c("Study.series", "Brain.region", "Subregion", "Sex", "METTL3",	"METTL14",	"VIRMA",	
        "WTAP",	"ZCCHC4",	"METTL5",	"TRMT112",	
        "METTL16",	"RBM15", "ZC3H13",	"CBLL1",	"YTHDC1",	"YTHDC2",	"YTHDF1",	
        "YTHDF2",	"YTHDF3",	"FMR1",	"IGF2BP1",	"IGF2BP2",	"IGF2BP3",	"FTO",	"ALKBH5")
tmat<-tmat[,rseq]

brcol<-brewer.pal(n=12, name = "Paired")
brcol[13:15]<-brewer.pal(n=3, name = "Set3")
names(brcol)<-unique (sort (tmat$Brain.region))

sbrcol<-brewer.pal(n=8, name = "Set3")
#sbrcol[10]<-"black"
names(sbrcol)<-unique(sort(tmat$Subregion))
sbrcol["L5"]<-"#E7298A"

ha = HeatmapAnnotation(df = tmat[,2:4], col = list(Brain.region = brcol, Subregion = sbrcol,
                                                   Sex=c("F"="green", "M"="yellow", "Mix"="orange", "Unknown/NA"="black")))

h<-Heatmap(as.matrix(t(tmat[,5:ncol(tmat)])), cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = ha, 
           column_labels = tmat$Study, heatmap_legend_param = list (title="Log2FC", direction="horizontal", title_position="topcenter"),
           row_split = data.frame(c(rep("Writers", 11), rep("Readers", 9), rep("Erasers", 2))),
           column_split = data.frame(tmat$Brain.region), column_title = NULL)

tiff("Figure2.tif", res = 300, width = 2500, height = 2000)
draw(h, heatmap_legend_side="bottom")
dev.off()
















