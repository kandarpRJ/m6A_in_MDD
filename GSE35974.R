#GSE35974


library(GEOquery)
library(ggplot2)
library(sva)
library(limma)

setwd("~/data/MDD_review/metaanalysis")


run_pca <- function(exp_mat) {
    pca <- prcomp(x = t(exp_mat), retx = TRUE, center = TRUE, scale. = TRUE)
    return(pca)
}

gseid<-"GSE35974"
Sys.setenv("VROOM_CONNECTION_SIZE"=1000000)
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
x<-x[ , which(apply(x, 2, var) != 0)]
exprs(gse)<-x

gse.pheno<-pData(gse)
gse.pheno<-data.frame(gse.pheno[gse.pheno$`disease status:ch1`=="unaffected" | gse.pheno$`disease status:ch1`=="depression",c("geo_accession", "characteristics_ch1.1", "characteristics_ch1",
                                                                                                                              "characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.4", 
                                                                                                                              "characteristics_ch1.5")])
colnames(gse.pheno)<-c("geo_accession", "label", "tissue", "age_at_death", "PMI", "sex", "brain_ph")
gse.pheno$tissue<-gsub("tissue: cerebellum", "CBLM", gse.pheno$tissue)
gse.pheno$label<-gsub("disease status: unaffected", "control", gse.pheno$label)
gse.pheno$label<-gsub("disease status: depression", "depression", gse.pheno$label)
gse.pheno$age<-gsub("age: ", "", gse.pheno$age)
gse.pheno$PMI<-gsub ("post-mortem interval \\(pmi\\): ", "", gse.pheno$PMI)
gse.pheno$sex<-gsub ("Sex: ", "", gse.pheno$sex)
gse.pheno$brain_ph<-gsub("ph: ", "", gse.pheno$brain_ph)
gse<-gse[,gse.pheno$geo_accession]


gse.pheno<-gse.pheno[pca_raw$x[,2] > -190,]


mod = model.matrix(~0+label+sex+age, data=gse.pheno)
mod0 = model.matrix(~1,data=gse.pheno)
svobj = sva(exprs(gse),mod,mod0)

#glm.sv1 <- glm(svobj$sv[,1]~gse.pheno[,"batch"])
#summary(glm.sv1)

gse.raw<-exprs(gse)
gse.sva<-removeBatchEffect(x = gse, covariates = svobj$sv, design = mod)
# 
exprs(gse)<-gse.sva
dat_lm<-foreach (i = 1:nrow(gse.raw), .combine=rbind) %dopar% {
    lm_model <- lm(gse.raw[i, ] ~ label+sex+age+PMI+brain_ph, data = gse.pheno)
    resids <- residuals(lm_model)
    resids + lm_model$coefficients[1] # Add back intercept term
}

pca_raw <- run_pca(gse.raw)
pca_lm <- run_pca(dat_lm)
pca_sva <- run_pca(gse.sva)


par(mfrow=c(1,2))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=gse.pheno$label))+geom_point()
ggplot(as.data.frame(pca_lm$x), aes(PC1, PC2, col=gse.pheno$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=gse.pheno$label))+geom_point()

design <- model.matrix(~0 + label + sex + age, gse.pheno)

fit <- lmFit(gse, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste("labeldepression", "labelcontrol", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

write.table(tT, file=paste(gseid,".txt",sep=""), row.names=F, sep="\t")








