rm(list = ls())
options(stringsAsFactors = FALSE)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(ggsignif)
library(limma)
library(xlsx)
library(ggpubr)
library(DESeq2)
load("NCBI_ID_sub.RData")
lncRNA = read.table("lncRNA.txt")[,1]##	3569
ensemble.lnc = read.table("Homo_sapiens.GRCh38.105.lncRNA.txt")
load("lncipedia_use.RData")
load("lncRNA_GSEA_res.RData")
load("KEGG.RData")
pathway = toupper(rownames(KEGG))

p.thres = data.frame(statistic = c(rep("NOM p-val", 3),"FDR q-val", "FWER p-val"), thres = c(0.01, 0.05, 0.1, 0.25, 0.25))
p = 2
path = paste0("/xxx/p", p.thres[p, "thres"])
setwd(path)

##	lncRNADisease
setwd(paste(path, "lncRNADisease", sep = "/"))
load("lnc2Disease.final.RData")
load("lncRNADisease.disease.alias.RData")
disease.alias.df = data.frame()
for(i in names(disease.alias)){
disease.alias.df = rbind(disease.alias.df, data.frame(Pathway = i, Disease = disease.alias[[i]]))
}
disease.alias.df = na.omit(disease.alias.df)
lnc2DiseasePathway.final = merge(lnc2Disease.final, disease.alias.df, by = "Pathway")
save(lnc2DiseasePathway.final, file = "lnc2DiseasePathway.final.RData")
load("lnc2DiseasePathway.final.RData")
data_use = unique(lnc2DiseasePathway.final[, 1:2])

##	lnc2Cancer
setwd(paste(path, "lnc2Cancer", sep = "/"))
load("lnc2Cancer.final.RData")
load("lnc2Cancer.disease.alias.RData")
disease.alias.df = data.frame()
for(i in names(disease.alias)){
disease.alias.df = rbind(disease.alias.df, data.frame(Pathway = i, Disease = disease.alias[[i]]))
}
disease.alias.df = na.omit(disease.alias.df)
lnc2CancerPathway.final = merge(lnc2Cancer.final, disease.alias.df, by = "Pathway")
save(lnc2CancerPathway.final, file = "lnc2CancerPathway.final.RData")
load("lnc2CancerPathway.final.RData")
data_use = unique(lnc2CancerPathway.final[, 1:2])

##	lncRNADisease, lnc2Cancer
setwd(paste(path, "lncRNADisease_lnc2Cancer", sep = "/"))
LD_L2C = unique(rbind(lnc2DiseasePathway.final, lnc2CancerPathway.final))
save(LD_L2C, file = "LD_L2C.RData")
load("LD_L2C.RData")
data_use = unique(LD_L2C[, 1:2])

lncRNA_all = union(LD_L2C$ncRNA.Ensembl, lncRNA)


##	lncRNAWiki
LncRNAWiki = read.csv("LncRNAWiki_Function.csv",sep=",", header = TRUE)
LncRNAWiki_lnc = data.frame()
for(i in unique(LncRNAWiki$symbol)){
LncRNAWiki_lnc = rbind(LncRNAWiki_lnc, NCBI_ID_sub[NCBI_ID_sub$Symbol == i, -1])
}
LncRNAWiki_lnc = unique(na.omit(LncRNAWiki_lnc))
colnames(LncRNAWiki_lnc) = "symbol"

LncRNAWiki_d = unique(na.omit(LncRNAWiki[LncRNAWiki$biological_context == "Disease", c("symbol", "context_detail")]))
LncRNAWiki_d$context_detail = toupper(LncRNAWiki_d$context_detail)
LncRNAWiki_d = merge(LncRNAWiki_d, LncRNAWiki_lnc, by.x = "symbol")
colnames(LncRNAWiki_d) = c("ncRNA.Symbol", "disease", "ncRNA.Ensembl")
save(LncRNAWiki_d, file = "LncRNAWiki_disease.RData")
load("LncRNAWiki_disease.RData")

lncRNA.stat = read.xlsx("ncRNA_disease_stat_summary.xlsx", sheetIndex = 3, startRow = 2)
disease.alias = list()
i = 28
toupper(lncRNA.stat[i, "Pathway"])
keyword = "MALARI"
d = unique(LncRNAWiki_d[grepl(toupper(keyword), LncRNAWiki_d$disease), "disease"])
d
disease.alias[[toupper(lncRNA.stat[i, "Pathway"])]] = NA
disease.alias[[toupper(lncRNA.stat[i, "Pathway"])]]
length(disease.alias)
unlist(lapply(disease.alias, length))
load("disease.alias.RData")
disease.alias.df = data.frame()
for(i in names(disease.alias)){
disease.alias.df = rbind(disease.alias.df, data.frame(Pathway = i, disease = disease.alias[[i]]))
}
disease.alias.df = na.omit(disease.alias.df)
LncRNAWiki_d = merge(LncRNAWiki_d, disease.alias.df, by = "disease")
load("LncRNAWiki_disease_pathway.RData")
data_use = unique(LncRNAWiki_d[, 3:4])

##	all databases
data_use = unique(rbind(lnc2DiseasePathway.final[, c("Pathway", "ncRNA.Ensembl")], rbind(lnc2CancerPathway.final[, c("Pathway", "ncRNA.Ensembl")], LncRNAWiki_d[, c("Pathway", "ncRNA.Ensembl")])))
save(data_use, file = "lncRNA2Disease_all.RData")

##	TCGA lncRNA fpkm
projects = c("BLCA", "BRCA", "COAD", "KIRC", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")
load(paste(path,"lncRNA.TCGA.RData", sep = "/"))
res = list()
for(project in projects){
print(project)
load(paste(path, project, "expr.fpkm.RData", sep = "/"))
print(dim(expr.fpkm))
meta = data.frame(cases = colnames(expr.fpkm), samples = strsplit2(colnames(expr.fpkm), "-")[, 4])
print(table(meta$samples))
meta[grep("01A", meta$samples), "Condition"] = "Tumor"
meta[grep("11A", meta$samples), "Condition"] = "Normal"
print(table(meta$Condition))
meta = na.omit(meta)
expr.fpkm.lnc = expr.fpkm[lncRNA.TCGA, meta$cases]
print(dim(expr.fpkm.lnc))
print(length(unique(strsplit2(rownames(expr.fpkm.lnc), "[.]")[, 1])))
expr.fpkm.lnc = aggregate(expr.fpkm.lnc, by = list(strsplit2(rownames(expr.fpkm.lnc), "[.]")[, 1]), FUN = mean)
rownames(expr.fpkm.lnc) = expr.fpkm.lnc$Group.1
expr.fpkm.lnc = expr.fpkm.lnc[, -1]
print(dim(expr.fpkm.lnc))
expr.fpkm.lnc = expr.fpkm.lnc[-which(rowSums(expr.fpkm.lnc) == 0), ]
print(dim(expr.fpkm.lnc))
group = factor(meta$Condition)
design = model.matrix(~0+group)
colnames(design)=levels(group)
rownames(design)=colnames(expr.fpkm.lnc)
contrast.matrix <- makeContrasts(Tumor-Normal, levels=design)
fit <- lmFit(expr.fpkm.lnc, design)
fit <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit)
tT = topTable(fit2,adjust='fdr',coef=1,number=Inf,p.value=1)
res[[project]] = tT
}
save(res, file = "TCGA_FPKM_limma_DElnc_res.RData")
load("TCGA_FPKM_limma_DElnc_res.RData")
projects = c("BLCA", "BRCA", "COAD", "KIRC", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")
all(projects == names(res))
TP = toupper(c("Bladder cancer", "Breast cancer", "Colorectal cancer", "Renal cell carcinoma", "Hepatocellular carcinoma", "Non-small cell lung cancer", "Non-small cell lung cancer", "Prostate cancer", "Gastric cancer", "Thyroid cancer", "Endometrial cancer"))
names(TP) = projects
data_use = data.frame()
for(i in projects){
data_use = rbind(data_use, data.frame(Pathway = TP[i], ncRNA.Ensembl = rownames(subset(res[[i]], abs(logFC)>1&adj.P.Val<0.05))))
}
save(data_use, file = "TCGA_FPKM_limma_DElnc_data_use.RData")
load("TCGA_FPKM_limma_DElnc_data_use.RData")

##	overlap
lnc2D.stat = data.frame(disease.pathway = unique(data_use$Pathway), row.names = unique(data_use$Pathway))
N = length(unique(ensemble.lnc[,1]))##	17741
for(i in unique(data_use$Pathway)){
lncRNA.disease = unique(data_use[data_use$Pathway == i, "ncRNA.Ensembl"])
lnc2D.stat[i, "lncRNA.disease"] = length(lncRNA.disease)
lncRNA.pathway = pathway_lncRNA_list[[i]]
lnc2D.stat[i, "lncRNA.pathway"] = length(lncRNA.pathway)
lnc2D.stat[i, "lncRNA.overlap"] = length(intersect(lncRNA.disease, lncRNA.pathway))
lnc2D.stat[i, "HGT.pvalue"] = phyper(lnc2D.stat[i, "lncRNA.overlap"]-1,lnc2D.stat[i, "lncRNA.disease"], N-lnc2D.stat[i, "lncRNA.disease"], lnc2D.stat[i, "lncRNA.pathway"], lower.tail=F)
}
lnc2D.stat[lnc2D.stat[, "HGT.pvalue"] > 0.05, "disease.pathway"]
lnc2D.stat$HGT.fdr = p.adjust(lnc2D.stat[, "HGT.pvalue"], method = "fdr")
lnc2D.stat[lnc2D.stat[, "HGT.fdr"] > 0.05, "disease.pathway"]
lnc2D.stat.order = lnc2D.stat[order(-lnc2D.stat$lncRNA.disease), ]
write.xlsx(lnc2D.stat.order, "lnc2D_pathway_lncRNA_overlap_HGT.xlsx")


##	PCC
projects = c("LIHC","STAD","BRCA","GBM","BLCA","PAAD", "PRAD", "COAD", "KIRC")
rownames(KEGG) = toupper(rownames(KEGG))
N = 19982##	GENCODE
for(project in projects){
print(project)
setwd(paste("/XXX", project, sep = "/"))
for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
load(paste0("lncRNA_CoExpression_mRNA_pcc", pcc.threshold, ".RData"))
lnc.m.HGT = list()
lnc.m.pathway = list()
for(lnc in intersect(lncRNA_all, names(lnc.m))){
if(length(lnc.m[[lnc]]) != 0){
HGT.res = data.frame()
for(j in 1:nrow(KEGG)){
HGT.res[j, "pathway"] = rownames(KEGG)[j]
pathway.gene = na.omit(as.numeric(KEGG[j, ]))
HGT.res[j, "pathway.gene"] = length(pathway.gene)
x1 = length(intersect(lnc.m[[lnc]], pathway.gene))
HGT.res[j, "overlap.gene"] = x1
x2 = length(lnc.m[[lnc]])
x3 = length(pathway.gene)
HGT.res[j, "p.value"] = phyper(x1-1,x2, N-x2, x3, lower.tail=F)
}
HGT.res.order = HGT.res[order(HGT.res$p.value), ]
HGT.res.order$fdr = p.adjust(HGT.res.order$p.value, method="fdr")
lnc.m.HGT[[lnc]] = HGT.res.order
lnc.m.pathway[[lnc]] = HGT.res.order[HGT.res.order$p.value < 0.05, "pathway"]
}
}
save(lnc.m.HGT, file = paste0("lncRNA_CoExpression_mRNA_pcc", pcc.threshold, "_HGT.RData"))
save(lnc.m.pathway, file = paste0("lncRNA_CoExpression_mRNA_pcc", pcc.threshold, "_pathwayList.RData"))
}
}

projects = c("BRCA", "KIRC", "BLCA", "COAD", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "STAD", "THCA", "UCEC")
TP = toupper(c("Breast cancer", "Renal cell carcinoma", "Bladder cancer", "Colorectal cancer", "Hepatocellular carcinoma", "Non-small cell lung cancer", "Non-small cell lung cancer", "Pancreatic cancer", "Prostate cancer", "Gastric cancer", "Thyroid cancer", "Endometrial cancer"))


##	overlap
projects = c("BLCA", "BRCA", "COAD", "GBM", "KIRC", "LAML", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "SKCM", "STAD", "THCA", "UCEC")
TP = toupper(c("Bladder cancer", "Breast cancer", "Colorectal cancer", "GLIOMA", "Renal cell carcinoma", "Acute Myeloid Leukemia", "Hepatocellular carcinoma", "Non-small cell lung cancer", "Non-small cell lung cancer", "Pancreatic cancer", "Prostate cancer", "MELANOMA", "Gastric cancer", "Thyroid cancer", "Endometrial cancer"))
names(TP) = projects
lnc2D.stat.order = xlsx::read.xlsx("lnc2D_pathway_lncRNA_overlap_HGT.xlsx", sheetIndex = 1)
lnc2D.stat.order = read.xlsx("lnc2D_pathway_lncRNA_overlap_HGT.xlsx")
lnc2D.stat.order = subset(lnc2D.stat.order, lncRNA.disease>20)
res.p = data.frame(ncFN = lnc2D.stat.order$HGT.pvalue, row.names = lnc2D.stat.order$disease.pathway)
res.fdr = data.frame(ncFN = lnc2D.stat.order$HGT.fdr, row.names = lnc2D.stat.order$disease.pathway)
N = length(unique(ensemble.lnc[,1]))##	17741

for(dis in c("100Kb", "500Kb", "1Mb")){
print(dis)
load(paste("/XXX/lncRNA_Colocation_pathway", dis, "mRNA.RData", sep = "_"))
lnc2D.stat = data.frame(disease.pathway = rownames(res.p), row.names = rownames(res.p))
for(i in rownames(res.p)){
lncRNA.disease = unique(data_use[data_use$Pathway == i, "ncRNA.Ensembl"])
lnc2D.stat[i, "lncRNA.disease"] = length(lncRNA.disease)
lncRNA.pathway = lnc.neighbor.m.pathway.list[[i]]
lnc2D.stat[i, "lncRNA.pathway"] = length(lncRNA.pathway)
lnc2D.stat[i, "lncRNA.overlap"] = length(intersect(lncRNA.disease, lncRNA.pathway))
lnc2D.stat[i, "HGT.pvalue"] = phyper(lnc2D.stat[i, "lncRNA.overlap"]-1,lnc2D.stat[i, "lncRNA.disease"], N-lnc2D.stat[i, "lncRNA.disease"], lnc2D.stat[i, "lncRNA.pathway"], lower.tail=F)
}
lnc2D.stat$HGT.fdr = p.adjust(lnc2D.stat[, "HGT.pvalue"], method = "fdr")
res.p[, paste0("CoLocation.", dis)] = lnc2D.stat$HGT.pvalue
res.fdr[, paste0("CoLocation.", dis)] = lnc2D.stat$HGT.fdr
}

for(dis in c("100Kb", "500Kb", "1Mb")){
print(dis)
load(paste("lncRNA_Colocation_pathway", dis, "mRNA.RData", sep = "_"))
print(table(unlist(lapply(lnc.neighbor.m.pathway.list, length))))
}

##	db
for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
for(i in c(1:7, 10:15)){
load(paste0("/XXX/", names(TP)[i], "/lncRNA_CoExpression_mRNA_pcc", pcc.threshold, "_mRNAList.RData"))
lncRNA.disease = unique(data_use[data_use$Pathway == TP[i], "ncRNA.Ensembl"])
lncRNA.pathway = lnc.m.pathway.list[[TP[i]]]
res.p[TP[i], paste0("CoExpression.pcc", pcc.threshold)] = phyper(length(intersect(lncRNA.disease, lncRNA.pathway))-1,length(lncRNA.disease), N-length(lncRNA.disease), length(lncRNA.pathway), lower.tail=F)
}
}


for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
for(i in c(1:7, 10:15)){
load(paste0("/XXX/", names(TP)[TP==i], "/lncRNA_CoExpression_mRNA_pcc", pcc.threshold, "_mRNAList.RData"))
lncRNA.disease = unique(data_use[data_use$Pathway == i, "ncRNA.Ensembl"])
lncRNA.pathway = lnc.m.pathway.list[[i]]
res.p[i, paste0("CoExpression.pcc", pcc.threshold)] = phyper(length(intersect(lncRNA.disease, lncRNA.pathway))-1,length(lncRNA.disease), N-length(lncRNA.disease), length(lncRNA.pathway), lower.tail=F)
}
}

lncRNA.disease = unique(data_use[data_use$Pathway == "NON-SMALL CELL LUNG CANCER", "ncRNA.Ensembl"])
for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
load(paste0("/XXX/lncRNA_CoExpression_mRNA_pcc", pcc.threshold, "_mRNAList.RData"))
lnc.m.pathway.list.LUAD = lnc.m.pathway.list
load(paste0("/XXX/lncRNA_CoExpression_mRNA_pcc", pcc.threshold, "_mRNAList.RData"))
lnc.m.pathway.list.LUSC = lnc.m.pathway.list
lncRNA.pathway = union(lnc.m.pathway.list.LUAD[["NON-SMALL CELL LUNG CANCER"]], lnc.m.pathway.list.LUSC[["NON-SMALL CELL LUNG CANCER"]])
res.p["NON-SMALL CELL LUNG CANCER", paste0("CoExpression.pcc", pcc.threshold)] = phyper(length(intersect(lncRNA.disease, lncRNA.pathway))-1,length(lncRNA.disease), N-length(lncRNA.disease), length(lncRNA.pathway), lower.tail=F)
}
save(res.p, file = "lncRNA2Disease_all_overlap_p.RData")


##	rank
projects = c("BLCA", "BRCA", "COAD", "GBM", "KIRC", "LAML", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "SKCM", "STAD", "THCA", "UCEC")
TP = toupper(c("Bladder cancer", "Breast cancer", "Colorectal cancer", "GLIOMA", "Renal cell carcinoma", "Acute Myeloid Leukemia", "Hepatocellular carcinoma", "Non-small cell lung cancer", "Non-small cell lung cancer", "Pancreatic cancer", "Prostate cancer", "MELANOMA", "Gastric cancer", "Thyroid cancer", "Endometrial cancer"))
names(TP) = projects
res = data.frame()
for(project in projects){
print(project)
lncRNA_d = data_use[data_use$Pathway == TP[project], "ncRNA.Ensembl"]
print(length(lncRNA_d))
res = rbind(res, data.frame(project = project, method = "ncFN", lncRNA = intersect(lncRNA_d, names(lncRNA_GSEA_res)), TPRank = NA))
for(lnc in intersect(lncRNA_d, names(lncRNA_GSEA_res))){
res[which(res$project == project & res$method == "ncFN" & res$lncRNA == lnc), "TPRank"] = match(TP[project], lncRNA_GSEA_res[[lnc]]$NAME)
}
for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
load(paste0("/XXX/", project, "/lncRNA_CoExpression_mRNA_pcc", pcc.threshold, "_HGT.RData"))
res = rbind(res, data.frame(project = project, method = paste0("CoExpression.pcc", pcc.threshold), lncRNA = intersect(lncRNA_d, names(lnc.m.HGT)), TPRank = NA))
for(lnc in intersect(lncRNA_d, names(lnc.m.HGT))){
res[which(res$project == project & res$method == paste0("CoExpression.pcc", pcc.threshold) & res$lncRNA == lnc), "TPRank"] = match(TP[project], lnc.m.HGT[[lnc]]$pathway)
}
}
for(dis in c("100Kb", "500Kb", "1Mb")){
print(dis)
load(paste("/XXX/lncRNA_Colocation_mRNA", dis, "HGT.RData", sep = "_"))
res = rbind(res, data.frame(project = project, method = paste0("CoLocation.", dis), lncRNA = intersect(lncRNA_d, names(lnc.neighbor.m.HGT)), TPRank = NA))
for(lnc in intersect(lncRNA_d, names(lnc.neighbor.m.HGT))){
res[which(res$project == project & res$method == paste0("CoLocation.", dis) & res$lncRNA == lnc), "TPRank"] = match(TP[project], lnc.neighbor.m.HGT[[lnc]]$pathway)
}
}
}
res$method = factor(res$method, levels = unique(res$method))


head(res)
res = data.frame()
for(project in unique(data_use$Pathway)){
print(project)
lncRNA_d = data_use[data_use$Pathway == project, "ncRNA.Ensembl"]
res = rbind(res, data.frame(project = project, method = "ncFN", lncRNA = intersect(lncRNA_d, names(lncRNA_GSEA_res)), TPRank = NA))
for(lnc in intersect(lncRNA_d, names(lncRNA_GSEA_res))){
res[which(res$project == project & res$method == "ncFN" & res$lncRNA == lnc), "TPRank"] = match(project, lncRNA_GSEA_res[[lnc]]$NAME)
}
if(project %in% TP){
for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
load(paste0("/XXX/", names(TP)[match(project, TP)], "/lncRNA_CoExpression_mRNA_pcc", pcc.threshold, "_HGT.RData"))
res = rbind(res, data.frame(project = project, method = paste0("CoExpression.pcc", pcc.threshold), lncRNA = intersect(lncRNA_d, names(lnc.m.HGT)), TPRank = NA))
for(lnc in intersect(lncRNA_d, names(lnc.m.HGT))){
res[which(res$project == project & res$method == paste0("CoExpression.pcc", pcc.threshold) & res$lncRNA == lnc), "TPRank"] = match(project, lnc.m.HGT[[lnc]]$pathway)
}
}
}
for(dis in c("100Kb", "500Kb", "1Mb")){
print(dis)
load(paste("/XXX/lncRNA_Colocation_mRNA", dis, "HGT.RData", sep = "_"))
res = rbind(res, data.frame(project = project, method = paste0("CoLocation.", dis), lncRNA = intersect(lncRNA_d, names(lnc.neighbor.m.HGT)), TPRank = NA))
for(lnc in intersect(lncRNA_d, names(lnc.neighbor.m.HGT))){
res[which(res$project == project & res$method == paste0("CoLocation.", dis) & res$lncRNA == lnc), "TPRank"] = match(project, lnc.neighbor.m.HGT[[lnc]]$pathway)
}
}
}
res$method = factor(res$method, levels = unique(res$method))
save(res, file = "lncRNA2Disease_all_rank.RData")

