#GSE17440

library(GEOquery)
library(ggplot2)
library(sva)
library(limma)

setwd("~/data/MDD_review/metaanalysis")


run_pca <- function(exp_mat) {
    pca <- prcomp(x = t(exp_mat), retx = TRUE, center = TRUE, scale. = TRUE)
    return(pca)
}

Sys.setenv("VROOM_CONNECTION_SIZE"=1000000)

gseid<-"GSE17440"
gse <- getGEO(gseid, GSEMatrix =TRUE, AnnotGPL=FALSE)
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

gse.pheno<-pData(gse)
gse.pheno<-data.frame(gse.pheno[, c("geo_accession", "latest mdd visit:ch1", "gender:ch1", "age:ch1", "tissue:ch1")])
colnames(gse.pheno)<-c("geo_accession", "label", "sex", "age", "tissue")
gse.pheno$label[c(1,2,4,8)]<-"depression"
gse.pheno$label[c(3,5,6,7)]<-"control"
gse.pheno$tissue<-gsub("Brain, frontal cortex", "FC", gse.pheno$tissue)
gse.pheno<-gse.pheno[-c(6,8),]
gse<-gse[,rownames(gse.pheno)]

gse.raw<-exprs(gse)

mod = model.matrix(~as.factor(label), data=gse.pheno)
mod0 = model.matrix(~1,data=gse.pheno)
svobj = sva(exprs(gse),mod,mod0)

gse.sva<-removeBatchEffect(x = exprs(gse), covariates = svobj$sv, design = mod)
exprs(gse)<-gse.sva

pca_raw <- run_pca(gse.raw)
pca_sva <- run_pca(gse.sva)

par(mfrow=c(1,2))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=gse.pheno$label, label=gse.pheno$geo_accession))+geom_point()+ggrepel::geom_label_repel()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=gse.pheno$label))+geom_point()

design <- model.matrix(~label + 0, gse.pheno)

fit <- lmFit(gse, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste("labeldepression", "labelcontrol", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.Symbol","Gene.Title"))
tT<-tT[order(tT$P.Value, decreasing = F),]
colnames(tT)<-c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title")

write.table(tT, file=paste(gseid,".txt",sep=""), row.names=F, sep="\t")















