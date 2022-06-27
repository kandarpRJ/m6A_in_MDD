#GSE54565


library(GEOquery)
library(ggplot2)
library(sva)
library(limma)

setwd("~/data/MDD_review/metaanalysis")


run_pca <- function(exp_mat) {
    pca <- prcomp(x = t(exp_mat), retx = TRUE, center = TRUE, scale. = TRUE)
    return(pca)
}

gseid<-"GSE54565"
gse <- getGEO(gseid, GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gse) > 1) idx <- grep("GPL570", attr(gse, "names")) else idx <- 1
gse <- gse[[idx]]

# make proper column names to match toptable 
fvarLabels(gse) <- make.names(fvarLabels(gse))

ex <- exprs(gse)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gse) <- log2(ex) }

x<-exprs(gse)
x<-x[which(apply(x, 1, var)!=0),]
exprs(gse)<-x

gse.pheno<-pData(gse)
gse.pheno<-data.frame(gse.pheno[gse.pheno$`disease state:ch1`=="Control" | gse.pheno$`disease state:ch1`=="MDD case",c("geo_accession", "characteristics_ch1.1", "characteristics_ch1")])
colnames(gse.pheno)<-c("geo_accession", "label", "tissue")
gse.pheno$tissue<-gsub("tissue: brain anterior cingulate cortex", "ACC", gse.pheno$tissue)
gse.pheno$label<-gsub("disease state: Control", "control", gse.pheno$label)
gse.pheno$label<-gsub("disease state: MDD case", "depression", gse.pheno$label)
gse.pheno$batch<-"A"
gse.pheno[pca_raw$x[,1]<0,"batch"]<-"B"

gse<-gse[,gse.pheno$geo_accession]

mod = model.matrix(~as.factor(label)+batch, data=gse.pheno)
mod0 = model.matrix(~1,data=gse.pheno)
svobj = sva(x,mod,mod0)

gse.raw<-x
gse.sva<-removeBatchEffect(x = x, batch = gse.pheno$batch, covariates = svobj$sv, design = mod[,1:2])

pca_raw <- run_pca(gse.raw)
pca_sva <- run_pca(gse.sva)

par(mfrow=c(1,2))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=gse.pheno$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=gse.pheno$label))+geom_point()

design <- model.matrix(~label + batch + 0, gse.pheno)

fit <- lmFit(x, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste("labeldepression", "labelcontrol", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)

gpl570<-read.csv("~/Downloads/gpl570.txt", sep = "\t", row.names = 1, header = T)
colnames(gpl570)<-c("Gene.symbol","Gene.title")
tT<-merge(tT, gpl570, by=0, all.x = T)
colnames(tT)<-c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Gene.symbol", "Gene.title")

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
tT<-tT[order(tT$P.Value, decreasing = F),]


write.table(tT, file=paste(gseid,".txt",sep=""), row.names=F, sep="\t")








