#GSE101521


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

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl",host = "grch37.ensembl.org"))

run_pca <- function(exp_mat) {
    pca <- prcomp(x = t(exp_mat), retx = TRUE, center = TRUE, scale. = TRUE)
    return(pca)
}

gseid<-"GSE101521"
gse=getGEO(filename="GSE101521_series_matrix.csv")
gse.pheno<-pData(gse)

sra<-read.csv("GSE101521_SraRunTable.txt")
sra<-sra[,c("Run","Sample.Name")]
gse.pheno<-merge(gse.pheno, sra, by.x = "geo_accession", by.y="Sample.Name", all.x = T)

countmat<-read.csv("GSE101521_counts.csv", header = TRUE, row.names = 1)
metadata<-gse.pheno[,c("tissue:ch1", "diagnosis:ch1", "Sex:ch1", "age (yrs):ch1", "brain ph:ch1",
                       "rin:ch1", "pmi:ch1")]
rownames(metadata)<-gse.pheno$title
colnames(metadata)<-c("region","label","sex","age_at_death", "brain_ph", "rin", "PMI")
metadata$label<-gsub ("non-psychiatric controls \\(CON\\)", "control", metadata$label)
metadata$label<-gsub ("DSM-IV major depressive disorder suicides \\(MDD-S\\)", "depression", metadata$label)
metadata$label<-gsub ("DSM-IV major depressive disorder non-suicides \\(MDD\\)", "depression", metadata$label)
metadata$region<-gsub("dorsal lateral prefrontal cortex \\(Brodmann Area 9\\)", "DLPFC", metadata$region)

y$counts<-y$counts[,rownames(metadata)]
##### OFC

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
write.table(finalDE, paste(gseid, ".txt", ""), quote = F, sep = "\t", row.names = F)



















