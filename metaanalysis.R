setwd("~/data/MDD_review/metaanalysis/")

library(metaRNASeq)
library(ComplexHeatmap)
library(clusterProfiler)
library(RColorBrewer)
library(org.Hs.eg.db)

meta<-read.csv("MDD_study_metadata.csv")

sel<-"DLPFC"

files<-list.files(pattern = "GSE.*.txt")
files<-files[grep(paste0(meta$file[meta$Brain.region  %in% sel],collapse = "|"), files)]

myfiles<-lapply(files, read.delim, header = T, sep="\t", fill=T, na.strings = c("","NA"))
files<-gsub(".txt", "", files)
files<-gsub("_DLPFC", "", files)

myfiles<-lapply(myfiles, function(x) na.omit(x))

mycols<-lapply(myfiles, function (x) data.frame(uname=make.names(x$Gene.symbol, unique = T), Gene.symbol=x$Gene.symbol, P.Value=x$P.Value, 
                                                logFC=x$logFC))
allmat<-Reduce(function (x, y) merge (x, y, by=c("uname")), mycols)
head (allmat)

pmat<-allmat[,c(1, grep("P.Value", colnames(allmat)))]
fcmat<-allmat[,c(1, grep("logFC", colnames(allmat)))]

plist<-lapply(2:10, function(x) pmat[,x])
names(plist)<-files

fclist<-lapply(2:10, function(x) fcmat[,x])
names(fclist)<-files

DE <- mapply(plist, FUN=function(x) ifelse(x <= 0.05, 1, 0))
colnames(DE)<-files
rownames(DE)<-pmat$uname

signsFC <- mapply(fclist, FUN=function(x) sign(x))
sumsigns <- apply(signsFC,1,sum)
commonsgnFC <- ifelse(abs(sumsigns)==dim(signsFC)[2], sign(sumsigns),0)


nreps<-c(59, 43, 28, 30, 26, 47, 76, 72, 251)

fishcomb <- fishercomb(plist, BHth = 0.05)
hist(fishcomb$rawpval, breaks=100, col="grey", main="Fisher method",
     xlab = "Raw p-values (meta-analysis)")

invnormcomb <- invnorm(plist,nrep=nreps, BHth = 0.05)
hist(invnormcomb$rawpval, breaks=100, col="grey",
     main="Inverse normal method",
     xlab = "Raw p-values (meta-analysis)")

DEresults <- data.frame(DE,
                        "DE.fishercomb"=ifelse(fishcomb$adjpval<=0.05,1,0),
                        "DE.invnorm"=ifelse(invnormcomb$adjpval<=0.05,1,0),
                        "adjP.fishercomb"=fishcomb$adjpval,
                        "adjP.invnormcomb"=invnormcomb$adjpval)
head(DEresults)

unionDE <- unique(c(fishcomb$DEindices,invnormcomb$DEindices))
FC.selecDE <- data.frame(DEresults[unionDE,],do.call(cbind,fclist)[unionDE,],
                         signFC=commonsgnFC[unionDE])
keepDE <- FC.selecDE[which(abs(FC.selecDE$signFC)==1),]
conflictDE <- FC.selecDE[which(FC.selecDE$signFC == 0),]
dim(FC.selecDE)


m6acounts<-read.table("temporal/human_FC_m6a_containing_genes_with_m6a_counts.txt", header = F)
#m6acounts<-read.table("human_crbm_m6a_containing_genes_with_m6a_counts.txt", header = F)
colnames(m6acounts)<-c("V1", "m6A_count")

gls<-bitr(m6acounts$V1, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop=FALSE)

m6acounts<-merge(m6acounts, gls, by.x="V1", by.y="ENSEMBL", all.x=T)

DEresults<-merge(DEresults, m6acounts, by.x = 0, by.y = "SYMBOL", all.x = T)
DEresults$m6A_count[is.na(DEresults$m6A_count)]<-0

selDEres<-DEresults[rowSums(DEresults[,2:10]==0)<9,]
selDEres<-selDEres[order(selDEres$adjP.fishercomb, selDEres$adjP.invnormcomb, decreasing = F),]
rownames(selDEres)<-selDEres$Row.names


write.table(selDEres, "metaanalysis_sign_genes.txt", quote = F, sep = "\t")

