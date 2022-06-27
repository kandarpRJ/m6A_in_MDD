#GSE54571


library(GEOquery)
library(doParallel)
library(ggplot2)
library(sva)
library(limma)

registerDoParallel(6)

setwd("~/data/MDD_review/metaanalysis")


run_pca <- function(exp_mat) {
    # Runs principal components analysis.
    # Input:
    #  exp_mat: data frame or matrix of expression values.
    #           Rows are genes and samples are columns.
    
    # Including some default settings just to be explicit
    pca <- prcomp(x = t(exp_mat), retx = TRUE, center = TRUE, scale. = TRUE)
    return(pca)
}

gseid<-"GSE54571"
gse <- getGEO(gseid, GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gse) > 1) idx <- grep("GPL570", attr(gse, "names")) else idx <- 1
gse <- gse[[idx]]

# make proper column names to match toptable 
fvarLabels(gse) <- make.names(fvarLabels(gse))

gse.pheno<-pData(gse)
gse.pheno<-data.frame(gse.pheno[gse.pheno$`disease state:ch1`=="Control" | gse.pheno$`disease state:ch1`=="MDD case",c("geo_accession", "characteristics_ch1.1", "characteristics_ch1")])
colnames(gse.pheno)<-c("geo_accession", "label", "tissue")
gse.pheno$tissue<-gsub("tissue: brain anterior cingulate cortex", "ACC", gse.pheno$tissue)
gse.pheno$label<-gsub("disease state: Control", "control", gse.pheno$label)
gse.pheno$label<-gsub("disease state: MDD case", "depression", gse.pheno$label)
gse<-gse[,gse.pheno$geo_accession]
#write.csv(gse.pheno, paste(gseid, ".pheno.csv", sep = ""), quote = F, row.names = F)
#write.csv(gse.expr, paste(gseid, ".expr.csv", sep = ""), quote = F)

mod = model.matrix(~as.factor(label), data=gse.pheno)
mod0 = model.matrix(~1,data=gse.pheno)
svobj = sva(exprs(gse),mod,mod0)

#new.pheno<-cbind(gse.pheno, svobj$sv)
#colnames(new.pheno)<-c("geo_accession","label","tissue","SV1","SV2","SV3","SV4","SV5","SV6")

#glm.sv1 <- glm(new.pheno[,"SV1"]~new.pheno[,"label"]) 
#summary(glm.sv1)

#glm.sv2 <- glm(new.pheno[,"SV2"]~new.pheno[,"label"]) 
#summary(glm.sv2)

#glm.sv3 <- glm(new.pheno[,"SV3"]~new.pheno[,"label"]) 
#summary(glm.sv3)

#glm.sv4 <- glm(new.pheno[,"SV4"]~new.pheno[,"label"]) 
#summary(glm.sv4)

#glm.sv5 <- glm(new.pheno[,"SV5"]~new.pheno[,"label"]) 
#summary(glm.sv5)

#glm.sv6 <- glm(new.pheno[,"SV6"]~new.pheno[,"label"]) 
#summary(glm.sv6)

gse.raw<-exprs(gse)
exprs(gse)<-removeBatchEffect(x = gse, covariates = svobj$sv, design = mod)

pca_raw <- run_pca(gse.raw)
pca_sva <- run_pca(exprs(gse))

par(mfrow=c(1,2))
ggplot(as.data.frame(pca_raw$x), aes(PC1, PC2, col=gse.pheno$label))+geom_point()
ggplot(as.data.frame(pca_sva$x), aes(PC1, PC2, col=gse.pheno$label))+geom_point()



#dat_lm<-foreach (i = 1:nrow(gse.expr), .combine=rbind) %dopar% {
#    lm_model <- lm(gse.expr[i, ] ~ label+SV1+SV2+SV3+SV4+SV5+SV6, data = new.pheno)
#    resids <- residuals(lm_model)
#    resids + lm_model$coefficients[1] # Add back intercept term
#}

design <- model.matrix(~label + 0, gse.pheno)

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








