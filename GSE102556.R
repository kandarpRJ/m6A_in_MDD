#GSE102556


library(edgeR)
library(GEOquery)
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

gseid<-"GSE102556"
gse=getGEO(filename="GSE102556_series_matrix.csv")
gse.pheno<-pData(gse)

sra<-read.csv("GSE102556_SraRunTable.csv")
sra<-sra[,c("Run","Sample.Name")]
gse.pheno<-merge(gse.pheno, sra, by.x = "geo_accession", by.y="Sample.Name", all.x = T)

countmat<-read.table("GSE102556_counts.csv", header = TRUE, row.names = 1)
s79<-read.table("SRR6443679_count.csv", header = T)
s79<-data.frame(row.names = s79[,1], SRR6443679=s79[,7])

s80<-read.table("SRR6443680_count.csv", header = T)
s80<-data.frame(row.names = s80[,1], SRR6443680=s80[,7])

countmat<-merge(countmat, s79, by=0)
rownames(countmat)<-countmat[,1]
countmat<-countmat[,-1]
countmat<-merge(countmat, s80, by=0)
rownames(countmat)<-countmat[,1]
countmat<-countmat[,-1]

metadata<-gse.pheno[,c("tissue:ch1", "phenotype:ch1", "gender:ch1", "age:ch1", "ph:ch1",
                        "rin:ch1", "pmi:ch1")]
rownames(metadata)<-gse.pheno$Run
colnames(metadata)<-c("region","label","sex","age_at_death", "brain_ph", "rin", "PMI")
metadata$label<-gsub ("MDD", "depression", metadata$label)
metadata$label<-gsub ("CTRL", "control", metadata$label)
metadata$region<-gsub("Orbitofrontal \\(OFC; BA11\\)", "OFC", metadata$region)
metadata$region<-gsub("Dorsolateral prefrontal cortex \\(dlPFC; BA8/9\\)", "DLPFC", metadata$region)
metadata$region<-gsub("Cingulate gyrus 25 \\(Cg25\\)", "CG25", metadata$region)
metadata$region<-gsub("Anterior Insula \\(aINS\\)", "AINS", metadata$region)
metadata$region<-gsub("Nucleus Accumbens \\(Nac\\)", "NACC", metadata$region)
metadata$region<-gsub("Subiculum \\(Sub\\)", "SUB", metadata$region)

#y$counts<-y$counts[,rownames(metadata)]
##### OFC

y<-DGEList(countmat[,metadata$region=="OFC" & (metadata$label=="control" | metadata$label=="depression")])
keep<-filterByExpr(y)
y<-y[keep,]
y<-calcNormFactors(y)
dat_cpm <- cpm(y, log = T)

ymeta<-metadata[metadata$region=="OFC" & (metadata$label=="control" | metadata$label=="depression"),]

dat_cpm<-subset(dat_cpm, select=-c(SRR5961809))
ymeta<-ymeta[rownames(ymeta) %NOTIN% c("SRR5961809"),]


mod = model.matrix(~0+label+sex+age_at_death, data=ymeta)
mod0 = model.matrix(~1,data=ymeta)
svobj = sva(dat_cpm,mod,mod0)

dat_lm<-foreach (i = 1:nrow(dat_cpm), .combine=rbind) %dopar% {
    lm_model <- lm(dat_cpm[i, ] ~ sex+age_at_death, data = ymeta)
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
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=ymeta$label, label=rownames(ymeta)))+geom_point()+ggrepel::geom_label_repel()
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
write.table(finalDE, paste(gseid, "_OFC.txt", ""), quote = F, sep = "\t", row.names = F)



###################### nAcc ###########################

y<-DGEList(countmat[,metadata$region=="NACC" & (metadata$label=="control" | metadata$label=="depression")])
keep<-filterByExpr(y)
y<-y[keep,]
y<-calcNormFactors(y)
dat_cpm <- cpm(y, log = T)

ymeta<-metadata[metadata$region=="NACC" & (metadata$label=="control" | metadata$label=="depression"),]

#dat_cpm<-subset(dat_cpm, select=-c(SRR5961990, SRR5962033, SRR5962032, SRR5962031, SRR5962030))
#ymeta<-ymeta[rownames(ymeta) %NOTIN% c("SRR5961990", "SRR5962033", "SRR5962032", "SRR5962031", "SRR5962030"),]

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

#dat_sva<-subset(dat_sva, select=-c(SRR5962033, SRR5962032, SRR5962031, SRR5962030))
#ymeta<-ymeta[rownames(ymeta) %NOTIN% c("SRR5962033", "SRR5962032", "SRR5962031", 'SRR5962030'),]

pca_raw <- run_pca(dat_cpm)
pca_lm <- run_pca(dat_lm)
pca_sva <- run_pca(dat_sva)

par(mfrow=c(1,3))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=ymeta$label, label=rownames(ymeta)))+geom_point()+ggrepel::geom_label_repel()
ggplot(as.data.frame(pca_lm$x), aes(PC1, PC2, col=ymeta$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=ymeta$label, label=rownames(ymeta)))+geom_point()+ggrepel::geom_label_repel()

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
write.table(finalDE, paste(gseid, "_nACC.txt", ""), quote = F, sep = "\t", row.names = F)




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

dat_sva<-removeBatchEffect(x = dat_cpm, covariates = svobj$sv, design = mod)

pca_raw <- run_pca(dat_cpm)
pca_lm <- run_pca(dat_lm)
pca_sva <- run_pca(dat_sva)

par(mfrow=c(1,3))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=ymeta$label, label = rownames(pca_raw$x)))+geom_point()+ggrepel::geom_text_repel()
dat_cpm<-subset(dat_cpm, select=-c(SRR5961878, SRR5961873, SRR5961883, SRR5961882, SRR5961875))
ymeta<-ymeta[colnames(dat_cpm),]
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
write.table(finalDE, paste(gseid, "_DLPFC.txt", ""), quote = F, sep = "\t", row.names = F)




###################### CG25 ###########################

y<-DGEList(countmat[,metadata$region=="CG25" & (metadata$label=="control" | metadata$label=="depression")])
keep<-filterByExpr(y)
y<-y[keep,]
y<-calcNormFactors(y)
dat_cpm <- cpm(y, log = T)

ymeta<-metadata[metadata$region=="CG25" & (metadata$label=="control" | metadata$label=="depression"),]
ymeta$batch<-"A"
ymeta$batch[rownames(ymeta) %in% rownames(pca_raw$x[pca_raw$x[,1]>0,])]<-"B"

#mod = model.matrix(~0+label+sex+age_at_death+batch, data=ymeta)
#mod0 = model.matrix(~1,data=ymeta)
#svobj = sva(dat_cpm,mod,mod0,n.sv = 4)

dat_lm<-foreach (i = 1:nrow(dat_cpm), .combine=rbind) %dopar% {
    lm_model <- lm(dat_cpm[i, ] ~ label+sex, data = ymeta)
    resids <- residuals(lm_model)
    resids + lm_model$coefficients[1] # Add back intercept term
}

glm.sv1 <- glm(svobj$sv[,1]~ymeta[,"sex"]) 
summary(glm.sv1)

#dat_sva<-removeBatchEffect(x = dat_cpm, covariates = svobj$sv, design = mod)

pca_raw <- run_pca(dat_cpm)
pca_lm <- run_pca(dat_lm)
pca_sva <- run_pca(dat_sva)

par(mfrow=c(1,3))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=ymeta$label))+geom_point()
ggplot(as.data.frame(pca_lm$x), aes(PC1, PC2, col=ymeta$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=ymeta$label))+geom_point()

design <- model.matrix(~0 + label + sex + age_at_death + batch, ymeta)

fit <- lmFit(dat_cpm, design)  # fit linear model

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
write.table(finalDE, paste(gseid, "_CG25.txt", ""), quote = F, sep = "\t", row.names = F)


###################### AINS ###########################

y<-DGEList(countmat[,metadata$region=="AINS" & (metadata$label=="control" | metadata$label=="depression")])
keep<-filterByExpr(y)
y<-y[keep,]
y<-calcNormFactors(y)
dat_cpm <- cpm(y, log = T)

ymeta<-metadata[metadata$region=="AINS" & (metadata$label=="control" | metadata$label=="depression"),]

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

dat_sva<-removeBatchEffect(x = dat_cpm, covariates = svobj$sv, design = mod)

pca_raw <- run_pca(dat_cpm)
pca_lm <- run_pca(dat_lm)
pca_sva <- run_pca(dat_sva)

par(mfrow=c(1,3))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=ymeta$label, label=rownames(ymeta)))+geom_point()+ggrepel::geom_label_repel()
dat_cpm<-subset(dat_cpm, select=-c(SRR5961961))
ymeta<-ymeta[colnames(dat_cpm),]
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
write.table(finalDE, paste(gseid, "_AINS.txt", ""), quote = F, sep = "\t", row.names = F)




###################### SUB ###########################

y<-DGEList(countmat[,metadata$region=="SUB" & (metadata$label=="control" | metadata$label=="depression")])
keep<-filterByExpr(y)
y<-y[keep,]
y<-calcNormFactors(y)
dat_cpm <- cpm(y, log = T)

ymeta<-metadata[metadata$region=="SUB" & (metadata$label=="control" | metadata$label=="depression"),]

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
write.table(finalDE, paste(gseid, "_SUB.txt", ""), quote = F, sep = "\t", row.names = F)



