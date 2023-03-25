//installing edgeR
install.packages(edgeR)

//load data into egdeR
library(edgeR)
setwd("~/Documents/IBDP MATERIAL /core guides /Extended Essay/Extended Essay Data")
load("EERD 2.37.27 AM.RData")
head(EERD)
class(EERD)

//Creating DGEList 
EERDGroups <- c("PSCLC", "PSCLC", "PSCLC", "PSCLC", "PSCLC", "PSCLC", "PSCLC", "PSCLC","PSCLC", "PSCLC", "PSCLC", "MSCLC", "MSLC", "MSCLC", "MSCLC", "MSCLC", "MSCLC", "MSCLC", "MSCLC", "MSCLC", "MSCLC", "MSCLC", "MSCLC", "MSCLC", "MSCLC", "MSCLC", "MSCLC")
d <- DEGList(counts = EERD, group = EERDGroups)
d

//Filtering the data 

head(d$counts)
head(cpm(d))
apply(d$counts, 2, sum)
keep <- rowSums(cpm(d)>=1)>=5  //specifying the limits to filter the CPM 
d <- d[keep,]
d$samples$lib.size <- colSums(d$counts)
d$samples

//Normalizing the data

di <- calcNormFactors(d)
di 

//Data Exploration

plotMDS(d, method = "bcv", col = as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col = 1:3, pch = 20)

//Estimating Dispersion 

//GLM estimates of dispersion 

design.mat <- model.matrix(~0+d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d, design.mat)
d2 <- estimateGLMTrendDisp(d2, design.mat, method = "power")
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2)

//Differential Gene Expression 

et12 <- exactTest(di, pair = c(1,2))
et13 <- exactTest(di, pair = c(2,1))
togTags(et12, n=10)
del <- row.names(di) [as.logical(del)]
plotSmear(et12, de.tags = deltags)
abline(h = c(-1,1), col = "black")
