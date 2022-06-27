#GSE53987

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

gseid<-"GSE53987"
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

gse.pheno<-pData(gse)
gse.pheno<-data.frame(gse.pheno[gse.pheno$`disease state:ch1`=="control" | gse.pheno$`disease state:ch1`=="major depressive disorder",
                                c("geo_accession", "disease state:ch1", "gender:ch1", "age:ch1", "tissue:ch1")])
colnames(gse.pheno)<-c("geo_accession", "label", "sex", "age_at_death", "tissue")
gse.pheno$tissue<-gsub("hippocampus", "HPC", gse.pheno$tissue)
gse.pheno$tissue<-gsub("Pre-frontal cortex \\(BA46\\)", "PFC", gse.pheno$tissue)
gse.pheno$tissue<-gsub("Associative striatum", "STR", gse.pheno$tissue)
gse.pheno$label<-gsub("major depressive disorder", "depression", gse.pheno$label)
gse<-gse[,rownames(gse.pheno)]

##############################################################################################################
## HPC


hpc<-gse[,gse.pheno$tissue=="HPC"]
hpc.pheno<-gse.pheno[gse.pheno$tissue=="HPC",]

hpc.raw<-exprs(hpc)

mod = model.matrix(~label + sex + age_at_death, data=hpc.pheno)
mod0 = model.matrix(~1,data=hpc.pheno)
svobj = sva(exprs(hpc),mod,mod0)

glm.sv1 <- glm(svobj$sv[,1]~hpc.pheno[,"sex"]) 
summary(glm.sv1)


hpc.sva<-removeBatchEffect(x = exprs(hpc), covariates = svobj$sv, design = mod)
exprs(hpc)<-hpc.sva

pca_raw <- run_pca(hpc.raw)
pca_sva <- run_pca(hpc.sva)

par(mfrow=c(1,2))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=hpc.pheno$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=hpc.pheno$label))+geom_point()

design <- model.matrix(~label + sex + age_at_death + 0, hpc.pheno)

fit <- lmFit(hpc, design)  # fit linear model

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

write.table(tT, file=paste(gseid,"_HPC.txt",sep=""), row.names=F, sep="\t")




##############################################################################################################
## PFC


pfc<-gse[,gse.pheno$tissue=="PFC"]
pfc.pheno<-gse.pheno[gse.pheno$tissue=="PFC",]

pfc.raw<-exprs(pfc)

mod = model.matrix(~label + sex + age_at_death, data=pfc.pheno)
mod0 = model.matrix(~1,data=pfc.pheno)
svobj = sva(exprs(pfc),mod,mod0)

glm.sv1 <- glm(svobj$sv[,1]~pfc.pheno[,"sex"]) 
summary(glm.sv1)


pfc.sva<-removeBatchEffect(x = exprs(pfc), covariates = svobj$sv, design = mod)
exprs(pfc)<-pfc.sva

pca_raw <- run_pca(pfc.raw)
pca_sva <- run_pca(pfc.sva)

par(mfrow=c(1,2))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=pfc.pheno$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=pfc.pheno$label))+geom_point()

design <- model.matrix(~label + sex + age_at_death + 0, pfc.pheno)

fit <- lmFit(pfc, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste("labeldepression", "labelcontrol", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
tT<-tT[order(tT$P.Value, decreasing = F),]
colnames(tT)<-c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title")

write.table(tT, file=paste(gseid,"_PFC.txt",sep=""), row.names=F, sep="\t")



##############################################################################################################
## STR


str<-gse[,gse.pheno$tissue=="STR"]
str.pheno<-gse.pheno[gse.pheno$tissue=="STR",]

str.raw<-exprs(str)

mod = model.matrix(~label + sex + age_at_death, data=str.pheno)
mod0 = model.matrix(~1,data=str.pheno)
svobj = sva(exprs(str),mod,mod0)

glm.sv1 <- glm(svobj$sv[,1]~str.pheno[,"sex"]) 
summary(glm.sv1)



str.sva<-removeBatchEffect(x = exprs(str), covariates = svobj$sv, design = mod)
exprs(str)<-str.sva

pca_raw <- run_pca(str.raw)
pca_sva <- run_pca(str.sva)

par(mfrow=c(1,2))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=str.pheno$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=str.pheno$label))+geom_point()

design <- model.matrix(~label + sex + age_at_death + 0, str.pheno)

fit <- lmFit(str, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste("labeldepression", "labelcontrol", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
tT<-tT[order(tT$P.Value, decreasing = F),]
colnames(tT)<-c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title")

write.table(tT, file=paste(gseid,"_STR.txt",sep=""), row.names=F, sep="\t")
















