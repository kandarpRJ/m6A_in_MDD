#GSE80655

library(GEOquery)
library(edgeR)
library(doParallel)
library(ggplot2)
library(sva)
library(limma)
library(tibble)
library(biomaRt)

setwd("~/data/MDD_review/metaanalysis")

registerDoParallel(6)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

run_pca <- function(exp_mat) {
    pca <- prcomp(x = t(exp_mat), retx = TRUE, center = TRUE, scale. = TRUE)
    return(pca)
}

gseid<-"GSE80655"
gse=getGEO(filename="GSE80655_series_matrix.csv")
gse.pheno<-pData(gse)

countmat<-read.table("GSE80655_counts.csv", header = TRUE, row.names = 1)
metadata<-gse.pheno[,c("brain region:ch1", "clinical diagnosis:ch1", "gender:ch1", "age at death:ch1", "brain ph:ch1",
                       "ethnicity:ch1", "post-mortem interval:ch1")]
rownames(metadata)<-gse.pheno$description
colnames(metadata)<-c("region","label","sex","age_at_death","brain_ph","ethnicity","PMI")
metadata$label<-gsub ("Major Depression", "depression", metadata$label)
metadata$label<-gsub ("Control", "control", metadata$label)


y<-DGEList(countmat[,metadata$region=="AnCg" & (metadata$label=="control" | metadata$label=="depression")])
keep<-filterByExpr(y)
y<-y[keep,]
y<-calcNormFactors(y)
dat_cpm <- cpm(y, log = T)

ymeta<-metadata[metadata$region=="AnCg" & (metadata$label=="control" | metadata$label=="depression"),]

mod = model.matrix(~0+label+sex+age_at_death, data=ymeta)
mod0 = model.matrix(~1,data=ymeta)
svobj = sva(dat_cpm,mod,mod0)

dat_lm<-foreach (i = 1:nrow(dat_cpm), .combine=rbind) %dopar% {
    lm_model <- lm(dat_cpm[i, ] ~ label+sex, data = ymeta)
    resids <- residuals(lm_model)
    resids + lm_model$coefficients[1] # Add back intercept term
}
#dat_lm<-dat_lm[ , which(apply(dat_lm, 2, var) != 0)]


glm.sv1 <- glm(svobj$sv[,1]~ymeta[,"sex"]) 
summary(glm.sv1)

dat_sva<-removeBatchEffect(x = dat_cpm, covariates = svobj$sv[,-c(1,6)], design = mod)

pca_raw <- run_pca(dat_cpm)
pca_lm <- run_pca(dat_lm)
pca_sva <- run_pca(dat_sva)

par(mfrow=c(1,3))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=ymeta$label))+geom_point()
ggplot(as.data.frame(pca_lm$x), aes(PC1, PC2, col=ymeta$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=ymeta$label))+geom_point()

design <- model.matrix(~0 + label + sex + age_at_death, ymeta)

fit <- lmFit(dat_sva, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste("labeldepression", "labelcontrol", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)

tT<-rownames_to_column(tT)

G_list <- getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
                values = tT$rowname, mart = mart)
finalDE<-merge(tT, G_list, by.x = "rowname", by.y = "ensembl_gene_id", all.x = T)
finalDE<-finalDE[,c(1,6,5,4,7,2,8,9)]
colnames(finalDE)<-c("ID", "adj.Pval", "P.Value", "t", "B", "logFC", "Gene.symbol", "Gene.title")
write.table(finalDE, "GSE80655_AnCg.txt", quote = F, sep = "\t", row.names = F)



###################### nAcc ###########################

y<-DGEList(countmat[,metadata$region=="nAcc" & (metadata$label=="control" | metadata$label=="depression")])
keep<-filterByExpr(y)
y<-y[keep,]
y<-calcNormFactors(y)
dat_cpm <- cpm(y, log = T)

ymeta<-metadata[metadata$region=="nAcc" & (metadata$label=="control" | metadata$label=="depression"),]

mod = model.matrix(~0+label+sex+age_at_death, data=ymeta)
mod0 = model.matrix(~1,data=ymeta)
svobj = sva(dat_cpm,mod,mod0)

dat_lm<-foreach (i = 1:nrow(dat_cpm), .combine=rbind) %dopar% {
    lm_model <- lm(dat_cpm[i, ] ~ label+sex, data = ymeta)
    resids <- residuals(lm_model)
    resids + lm_model$coefficients[1] # Add back intercept term
}
#dat_lm<-dat_lm[ , which(apply(dat_lm, 2, var) != 0)]


glm.sv1 <- glm(svobj$sv[,1]~ymeta[,"sex"]) 
summary(glm.sv1)

dat_sva<-removeBatchEffect(x = dat_cpm, covariates = svobj$sv, design = mod)

pca_raw <- run_pca(dat_cpm)
pca_lm <- run_pca(dat_lm)
pca_sva <- run_pca(dat_sva)

par(mfrow=c(1,3))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=ymeta$label))+geom_point()
ggplot(as.data.frame(pca_lm$x), aes(PC1, PC2, col=ymeta$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=ymeta$label))+geom_point()

design <- model.matrix(~0 + label + sex + age_at_death, ymeta)

fit <- lmFit(dat_sva, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste("labeldepression", "labelcontrol", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)

tT<-rownames_to_column(tT)

G_list <- getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
                values = tT$rowname, mart = mart)
finalDE<-merge(tT, G_list, by.x = "rowname", by.y = "ensembl_gene_id", all.x = T)
finalDE<-finalDE[,c(1,6,5,4,7,2,8,9)]
colnames(finalDE)<-c("ID", "adj.Pval", "P.Value", "t", "B", "logFC", "Gene.symbol", "Gene.title")
finalDE<-finalDE[order(finalDE$adj.Pval, decreasing = F),]
write.table(finalDE, "GSE80655_nACC.txt", quote = F, sep = "\t", row.names = F)




###################### DLPFC ###########################

y<-DGEList(countmat[,metadata$region=="DLPFC" & (metadata$label=="control" | metadata$label=="depression")])
keep<-filterByExpr(y)
y<-y[keep,]
y<-calcNormFactors(y)
dat_cpm <- cpm(y, log = T)

ymeta<-metadata[metadata$region=="DLPFC" & (metadata$label=="control" | metadata$label=="depression"),]

mod = model.matrix(~0+label+sex+age_at_death, data=ymeta)
mod0 = model.matrix(~1,data=ymeta)
svobj = sva(dat_cpm,mod,mod0)

dat_lm<-foreach (i = 1:nrow(dat_cpm), .combine=rbind) %dopar% {
    lm_model <- lm(dat_cpm[i, ] ~ label+sex, data = ymeta)
    resids <- residuals(lm_model)
    resids + lm_model$coefficients[1] # Add back intercept term
}

glm.sv1 <- glm(svobj$sv[,1]~ymeta[,"sex"]) 
summary(glm.sv1)

dat_sva<-removeBatchEffect(x = dat_cpm, covariates = svobj$sv[,-2], design = mod)

pca_raw <- run_pca(dat_cpm)
pca_lm <- run_pca(dat_lm)
pca_sva <- run_pca(dat_sva)

par(mfrow=c(1,3))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=ymeta$label))+geom_point()
ggplot(as.data.frame(pca_lm$x), aes(PC1, PC2, col=ymeta$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=ymeta$label))+geom_point()

design <- model.matrix(~0 + label + sex + age_at_death, ymeta)

fit <- lmFit(dat_sva, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste("labeldepression", "labelcontrol", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)

tT<-rownames_to_column(tT)

G_list <- getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
                values = tT$rowname, mart = mart)
finalDE<-merge(tT, G_list, by.x = "rowname", by.y = "ensembl_gene_id", all.x = T)
finalDE<-finalDE[,c(1,6,5,4,7,2,8,9)]
colnames(finalDE)<-c("ID", "adj.Pval", "P.Value", "t", "B", "logFC", "Gene.symbol", "Gene.title")
finalDE<-finalDE[order(finalDE$adj.Pval, decreasing = F),]
write.table(finalDE, "GSE80655_DLPFC.txt", quote = F, sep = "\t", row.names = F)






