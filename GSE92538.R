#GSE92538

library(GEOquery)
library(ggplot2)
library(sva)
library(limma)

setwd("~/data/MDD_review/metaanalysis")


run_pca <- function(exp_mat) {
    pca <- prcomp(x = t(exp_mat), retx = TRUE, center = TRUE, scale. = TRUE)
    return(pca)
}

Sys.setenv("VROOM_CONNECTION_SIZE"=10000000)

gseid<-"GSE92538"
gse <- getGEO(gseid, GSEMatrix =TRUE, AnnotGPL=TRUE)
#if (length(gse) > 1) idx <- grep("GPL570", attr(gse, "names")) else idx <- 1
gse <- gse[[1]]

#make proper column names to match toptable
fvarLabels(gse) <- make.names(fvarLabels(gse))

ex <- exprs(gse)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gse) <- log2(ex) }
gpl<-getGEO("GPL10526")

dat_x<-exprs(gse)
dat_x<-merge(dat_x, Table(gpl)[,1:2], by.x=0, by.y="ID", all.x=T)
dat_1<-Reduce(function (x, y) merge(x, y, by="SYMBOL"), lapply(dat_x[,2:129], function (x) aggregate(x~SYMBOL, dat_x, mean)))
colnames(dat_1)<-c("Gene.symbol", colnames(dat_x[,2:129]))
dat_1<-sapply(dat_1[,2:ncol(dat_1)], function(x) (x-mean(x))/sd(x))
rownames(dat_1)<-sort (unique (dat_x$Gene.symbol))

gse.pheno<-pData(gse)

gse.pheno<-data.frame(gse.pheno[gse.pheno$`diagnosis:ch1`=="Control" | gse.pheno$`diagnosis:ch1`=="Major Depressive Disorder",
                                c("geo_accession", "diagnosis:ch1", "gender:ch1", "age:ch1", "tissue:ch1")])
colnames(gse.pheno)<-c("geo_accession", "label", "sex", "age_at_death", "tissue")
gse.pheno$tissue<-gsub("Dorsolateral Prefrontal Cortex", "DLPFC", gse.pheno$tissue)
gse.pheno$label<-gsub("Major Depressive Disorder", "depression", gse.pheno$label)
gse.pheno$label<-gsub("Control", "control", gse.pheno$label)
gse.pheno$ArrayVersion<-"GPL10526"

gse <- getGEO(gseid, GSEMatrix =TRUE, AnnotGPL=TRUE)
# if (length(gse) > 1) idx <- grep("GPL570", attr(gse, "names")) else idx <- 1
gse <- gse[[2]]

#make proper column names to match toptable
fvarLabels(gse) <- make.names(fvarLabels(gse))

ex <- exprs(gse)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gse) <- log2(ex) }
#gpl<-getGEO("GPL17027")   no gene sym. using id to gene map of prev platform

dat_x<-exprs(gse)
dat_x<-merge(dat_x, Table(gpl)[,1:2], by.x=0, by.y="ID", all.x=T)
dat_2<-Reduce(function (x, y) merge(x, y, by="SYMBOL"), lapply(dat_x[,2:236], function (x) aggregate(x~SYMBOL, dat_x, mean)))
colnames(dat_2)<-c("Gene.symbol", colnames(dat_x[,2:236]))
dat_2<-sapply(dat_2[,2:ncol(dat_2)], function(x) (x-mean(x))/sd(x))
rownames(dat_2)<-sort (unique (dat_x$SYMBOL))





pheno<-pData(gse)

pheno<-data.frame(pheno[pheno$`diagnosis:ch1`=="Control" | pheno$`diagnosis:ch1`=="Major Depressive Disorder",
                                c("geo_accession", "diagnosis:ch1", "gender:ch1", "age:ch1", "tissue:ch1")])
colnames(pheno)<-c("geo_accession", "label", "sex", "age_at_death", "tissue")
pheno$tissue<-gsub("Dorsolateral Prefrontal Cortex", "DLPFC", pheno$tissue)
pheno$label<-gsub("Major Depressive Disorder", "depression", pheno$label)
pheno$label<-gsub("Control", "control", pheno$label)
pheno$ArrayVersion<-"GPL17027"
gse.pheno<-rbind(gse.pheno, pheno)

dat_mat<-cbind(dat_1, dat_2)
dat_mat<-dat_mat[,rownames(gse.pheno)]

mod = model.matrix(~as.factor(label)+sex+age_at_death+ArrayVersion, data=gse.pheno)
mod0 = model.matrix(~1,data=gse.pheno)
svobj = sva(dat_mat,mod,mod0)

gse.sva<-removeBatchEffect(x = dat_mat, batch = gse.pheno$ArrayVersion, covariates = svobj$sv, design = mod[,1:52])

pca_raw <- run_pca(dat_mat)
pca_sva <- run_pca(gse.sva)

par(mfrow=c(1,2))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=gse.pheno$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=gse.pheno$label))+geom_point()

design <- model.matrix(~label + sex + age_at_death + ArrayVersion + 0, gse.pheno)

fit <- lmFit(gse.sva, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste("labeldepression", "labelcontrol", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)
tT$ID<-rownames(tT)
tT$Gene.symbol<-rownames(tT)
tT$Gene.title<-rownames(tT)

tT<-tT[,c("ID", "adj.P.Val", "P.Value", "t", "B","logFC", "Gene.symbol", "Gene.title")]

tT<-tT[order(tT$P.Value, decreasing = F),]


write.table(tT, file=paste(gseid,".txt",sep=""), row.names=F, sep="\t")












