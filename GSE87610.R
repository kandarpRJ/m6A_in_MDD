#GSE87610

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

gseid<-"GSE87610"
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
gse.pheno<-data.frame(gse.pheno[gse.pheno$`genotype:ch1`=="Unaffected Comparison subject" | gse.pheno$`genotype:ch1`=="Major Depressive Disorder",
                                c("geo_accession", "genotype:ch1", "tissue:ch1", "description")])
colnames(gse.pheno)<-c("geo_accession", "label", "tissue", "description")
gse.pheno$tissue<-gsub("dorsolateral prefrontal cortex", "DLPFC", gse.pheno$tissue)
gse.pheno$label<-gsub("Unaffected Comparison subject", "control", gse.pheno$label)
gse.pheno$label<-gsub("Major Depressive Disorder", "depression", gse.pheno$label)
gse<-gse[,rownames(gse.pheno)]

################################################################################
## L3

l3.pheno<-gse.pheno[grep ("layer 3", gse.pheno$description),]
l3<-gse[,rownames(l3.pheno)]

l3.raw<-exprs(l3)

mod = model.matrix(~as.factor(label), data=l3.pheno)
mod0 = model.matrix(~1,data=l3.pheno)
svobj = sva(exprs(l3),mod,mod0)

l3.sva<-removeBatchEffect(x = exprs(l3), covariates = svobj$sv, design = mod)
exprs(l3)<-l3.sva

pca_raw <- run_pca(l3.raw)
pca_sva <- run_pca(l3.sva)

par(mfrow=c(1,2))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=l3.pheno$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=l3.pheno$label))+geom_point()

design <- model.matrix(~label + 0, l3.pheno)

fit <- lmFit(l3, design)  # fit linear model

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
tT<-na.omit(tT)
write.table(tT, file=paste(gseid,"_L3.txt",sep=""), row.names=F, sep="\t")


################################################################################
## L5

l5.pheno<-gse.pheno[grep ("layer 5", gse.pheno$description),]
l5<-gse[,rownames(l5.pheno)]

l5.raw<-exprs(l5)

mod = model.matrix(~as.factor(label), data=l5.pheno)
mod0 = model.matrix(~1,data=l5.pheno)
svobj = sva(exprs(l5),mod,mod0)

l5.sva<-removeBatchEffect(x = exprs(l5), covariates = svobj$sv, design = mod)
exprs(l5)<-l5.sva

pca_raw <- run_pca(l5.raw)
pca_sva <- run_pca(l5.sva)

par(mfrow=c(1,2))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=l5.pheno$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=l5.pheno$label))+geom_point()

design <- model.matrix(~label + 0, l5.pheno)

fit <- lmFit(l5, design)  # fit linear model

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
tT<-na.omit(tT)
write.table(tT, file=paste(gseid,"_L5.txt",sep=""), row.names=F, sep="\t")

































