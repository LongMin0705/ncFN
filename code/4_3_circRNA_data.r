rm(list = ls())
options(stringsAsFactors = FALSE)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(limma)
library(data.table)
library(RColorBrewer)
library(edgeR)

# lncRNADisease
lncRNADisease.circ = read.xlsx("experimental circRNA-disease information.xlsx")
lncRNADisease.circ.homo = lncRNADisease.circ[lncRNADisease.circ$Species == "Homo sapiens", ]
circ2Disease = unique(lncRNADisease.circ.homo[, c("ncRNA.Symbol", "Disease.Name")])
circ2Disease$Disease.Name = toupper(circ2Disease$Disease.Name)
load("lncRNADisease_circRNA_disease_alias.RData")
disease.alias.df = data.frame()
for(i in names(disease.alias)){
disease.alias.df = rbind(disease.alias.df, data.frame(Pathway = i, Disease = disease.alias[[i]]))
}
circ2Disease.dps = unique(circ2Disease[circ2Disease$Disease.Name %in% disease.alias.df$Disease, ])
circ2Disease.dps$Pathway = disease.alias.df[match(circ2Disease.dps$Disease.Name, disease.alias.df$Disease), "Pathway"]
circ2Disease.dps$ncRNA.id = circRNA_info[match(circ2Disease.dps$ncRNA.Symbol, circRNA_info$gene_symbol), "circRNA_ID"]
save(circ2Disease.dps, file = "circ2Disease.dps1.RData")
load("circ2Disease.dps1.RData")

idss = which(is.na(circ2Disease.dps$ncRNA.id))
unique(circ2Disease.dps[idss, "ncRNA.Symbol"])
circ2Disease.dps[, "ncRNA.Symbol"] = gsub("_circRNA_", "_circ_", circ2Disease.dps[, "ncRNA.Symbol"])
circ2Disease.dps[which(circ2Disease.dps$ncRNA.Symbol == "circRNA0003906"), "ncRNA.Symbol"] = "hsa_circ_0003906"
circ2Disease.dps[grep("^hsa_circ_", circ2Disease.dps$ncRNA.Symbol), "ncRNA.id"] = circ2Disease.dps[grep("^hsa_circ_", circ2Disease.dps$ncRNA.Symbol), "ncRNA.Symbol"]
circ2Disease.dps[which(circ2Disease.dps$ncRNA.id == "hsa_circ_0126897_CBC1"), "ncRNA.id"] = "hsa_circ_0126897"
circ2Disease.dps[which(nchar(circ2Disease.dps$ncRNA.id) != 16), "ncRNA.id"] = NA
save(circ2Disease.dps, file = "circ2Disease.dps2.RData")
load("circ2Disease.dps2.RData")

idss = which(is.na(circ2Disease.dps$ncRNA.id))
for(i in idss){
a = circ2Disease.dps[i, "ncRNA.Symbol"]
if(nchar(a) > 6){
b = substring(a, nchar(a)-5, nchar(a))
circ2Disease.dps[i, "ncRNA.Symbol"] = paste0("hsa_circ_", b)
}
}
circ2Disease.dps[idss, "ncRNA.Symbol"] %in% circID_to_name$circRNA_name
circ2Disease.dps[idss, "ncRNA.id"] = circID_to_name[match(circ2Disease.dps[idss, "ncRNA.Symbol"], circID_to_name$circRNA_name), "circRNA_ID"]
save(circ2Disease.dps, file = "circ2Disease.dps3.RData")
load("circ2Disease.dps3.RData")
idss = which(is.na(circ2Disease.dps$ncRNA.id))
unique(circ2Disease.dps[idss, "ncRNA.Symbol"])
circ2Disease.dps = unique(na.omit(circ2Disease.dps))
circ2Disease.final = unique(circ2Disease.dps[, 3:4])
save(circ2Disease.final, file = "circ2Disease.final.RData")
load("circ2Disease.final.RData")


##	circRNADisease
circRNADisease = xlsx::read.xlsx("circRNADisease-2017-12-25.xlsx", sheetIndex = 1)
circRNADisease$species = toupper(circRNADisease$species)
circRNADisease.sub = circRNADisease[circRNADisease$species=="HUMAN", c("circRNA.id", "circRNA.name", "circRNA.synonyms", "disease")]
idx = grep("hsa_circ_", circRNADisease.sub$circRNA.id)
circ2Disease = circRNADisease.sub[idx, c("circRNA.id", "disease")]
idx = grep("hsa_circ_", circRNADisease.sub$circRNA.name)
temp = circRNADisease.sub[idx, c("circRNA.name", "disease")]
colnames(temp) = colnames(circ2Disease)
circ2Disease = unique(rbind(circ2Disease, temp))
idx = grep("hsa_circ_", circRNADisease.sub$circRNA.synonyms)
temp = circRNADisease.sub[idx, c("circRNA.synonyms", "disease")]
colnames(temp) = colnames(circ2Disease)
circ2Disease = rbind(circ2Disease, temp)
circ2Disease = unique(circ2Disease)
circ2Disease$disease = toupper(circ2Disease$disease)
table(nchar(circ2Disease$circRNA.id))
length(intersect(circRNA, unique(circ2Disease$circRNA.id)))
circ2Disease[which(nchar(circ2Disease$circRNA.id) == 15), "circRNA.id"] %in%  circID_to_name$circRNA_name
circ2Disease[which(circ2Disease$circRNA.id %in% circID_to_name$circRNA_name), "circRNA.id"] = na.omit(circID_to_name[match(circ2Disease$circRNA.id, circID_to_name$circRNA_name), "circRNA_ID"])

circ2Disease.dps = circ2Disease[circ2Disease$circRNA.id %in% circRNA, ]
circ2Disease.dps[, "Pathway"] = circ2Disease.dps[, "disease"]
circ2Disease.dps[, "Pathway"] %in% pathway
circ2Disease.dps[8, "Pathway"] = NA
circ2Disease.dps[13, "Pathway"] = "THYROID CANCER"
circ2Disease.dps[18, "Pathway"] = NA
circ2Disease.dps[19, "Pathway"] = "LONG-TERM DEPRESSION"
circ2Disease.dps[21, "Pathway"] = NA
circ2Disease.dps[22, "Pathway"] = NA
circ2Disease.dps[25, "Pathway"] = "PANCREATIC CANCER"
circ2Disease.dps[26, "Pathway"] = "PANCREATIC CANCER"
circ2Disease.dps = unique(na.omit(circ2Disease.dps))
save(circ2Disease.dps, file = "circ2Disease.dps.RData")
load("circ2Disease.dps.RData")


##	circFunBase
circFB = read.table("circFunBase_Homo_sapiens_circ.txt", header = TRUE, sep = "\t", quote = "")
circFB$Function =toupper(circFB$Function)
length(grep("^hsa_circ_", circFB$circRNA))
for(i in setdiff(1:nrow(circFB), grep("^hsa_circ_", circFB$circRNA))){
if(!is.na(match(circFB[i, "Position"], circBank$position)))
circFB[i, "circRNA"] = circBank[match(circFB[i, "Position"], circBank$position), "circbaseID"]
}
length(grep("^hsa_circ_", circFB$circRNA))
circFB_sub = cbind(circFB[, "circRNA"], strsplit2(circFB$Function, "; "))
circFB_melt = data.frame()
for(i in 1:nrow(circFB_sub)){
circFB_melt = rbind(circFB_melt, data.frame(circRNA = circFB_sub[i, 1], disease = circFB_sub[i, -1]))
}
length(which(circFB_melt$disease == ""))
circFB_melt = circFB_melt[-which(circFB_melt$disease == ""), ]
circFB_use = unique(cbind(circFB_melt[, "circRNA"], disease = toupper(strsplit2(circFB_melt$disease, " [[]")[,1])))
circFB_use = data.frame(circRNA = circFB_use[, 1], disease = circFB_use[, 2])
length(intersect(circRNA, unique(circFB_use$circRNA)))

circRNA.stat = read.xlsx("ncRNA_disease_stat_summary.xlsx", sheet = "miRNA.HMDD", startRow = 2)
disease.alias = list()
i = 28
toupper(circRNA.stat[i, "Pathway"])
keyword = "MALARIA"
d = unique(circFB_use[grepl(toupper(keyword), circFB_use$disease), "disease"])
d
[-c(2,6,11,17,21,24,26,27,29)]
disease.alias[[toupper(circRNA.stat[i, "Pathway"])]] = NA
disease.alias[[toupper(circRNA.stat[i, "Pathway"])]]
length(disease.alias)
unlist(lapply(disease.alias, length))
load("circFunBase_disease.alias.RData")

disease.alias.df = data.frame()
for(i in names(disease.alias)){
disease.alias.df = rbind(disease.alias.df, data.frame(Pathway = i, Disease = disease.alias[[i]]))
}
disease.alias.df = na.omit(disease.alias.df)
circFB_dps = circFB_use[circFB_use$disease %in% disease.alias.df$Disease, ]
circFB_dps$Pathway = disease.alias.df[match(circFB_dps$disease, disease.alias.df$Disease), "Pathway"]
circFB_dps$circRNA = tolower(circFB_dps$circRNA)
length(grep("^hsa_circ_", circFB_dps$circRNA))
circFB_dps = circFB_dps[grep("^hsa_circ_", circFB_dps$circRNA), ]
circFB_dps = circFB_dps[nchar(circFB_dps$circRNA) == 16, ]
save(circFB_dps, file = "circFB_dps.RData")
load("circFB_dps.RData")


##	circad
setwd(path)
circad = read.table("circad.txt", sep = "\t", header = TRUE)
circad$circRNA.Name = tolower(circad$circRNA.Name)
length(grep("hsa_circ_", circad$circRNA.Name))
circad$Alias = tolower(circad$Alias)
length(grep("hsa_circ_", circad$Alias))
circad[grep("hsa_circ_", circad$circRNA.Name), "circRNA_ID"] = circad[grep("hsa_circ_", circad$circRNA.Name), "circRNA.Name"]
circad[is.na(circad$circRNA_ID) & grepl("^hsa_circ_", circad$Alias), "circRNA_ID"] = circad[is.na(circad$circRNA_ID) & grepl("^hsa_circ_", circad$Alias), "Alias"]
circad$Genome.Locus.new = paste(strsplit2(circad$Genome.Locus, ":")[, 1], strsplit2(circad$Genome.Locus, ":")[, 2], sep = ":")
circRNA_ID$Genome.Locus = paste(circRNA_ID$chrom, paste(circRNA_ID$start, circRNA_ID$end, sep = "-"), sep = ":")
length(which(circad$Genome.Locus.new %in% circRNA_ID$Genome.Locus))
circad$circRNA_ID = circRNA_ID[match(circad$Genome.Locus, circRNA_ID$Genome.Locus), "circRNA.ID"]
circad[is.na(circad$circRNA_ID), "circRNA_ID"] = circRNA_ID[match(circad[is.na(circad$circRNA_ID), "Genome.Locus.new"], circRNA_ID$Genome.Locus), "circRNA.ID"]
save(circad, file = "circad_transfer.RData")

circad_sub = unique(na.omit(circad[, c("circRNA_ID", "Disease")]))
circad_sub$Disease = toupper(circad_sub$Disease)
circad_sub = circad_sub[-94,]

circRNA.stat = read.xlsx("ncRNA_disease_stat_summary.xlsx", sheetIndex = 3, startRow = 2)
disease.alias = list()
i = 28
toupper(circRNA.stat[i, "Pathway"])
keyword = "MALARIA"
d = unique(circad_sub[grepl(toupper(keyword), circad_sub$Disease), "Disease"])
d
disease.alias[[toupper(circRNA.stat[i, "Pathway"])]] = NA
disease.alias[[toupper(circRNA.stat[i, "Pathway"])]]
length(disease.alias)
unlist(lapply(disease.alias, length))
save(disease.alias, file = "circad_disease.alias.RData")
disease.alias.df = data.frame()
for(i in names(disease.alias)){
disease.alias.df = rbind(disease.alias.df, data.frame(Pathway = i, Disease = disease.alias[[i]]))
}
disease.alias.df = na.omit(disease.alias.df)
circad_use = circad_sub[circad_sub$Disease %in% disease.alias.df$Disease, ]
circad_use$Pathway = disease.alias.df[match(circad_use$Disease, disease.alias.df$Disease), "Pathway"]circRNA, 32 disease, 17 disease pathway
circad_use = circad_use[grep("^hsa_circ_", circad_use$circRNA_ID), ]
circad_use = circad_use[nchar(circad_use$circRNA_ID) == 16, ]
save(circad_use, file = "circad_use.RData")
load("circad_use.RData")


##	all databases
load("circ2Disease.final.RData")##	lncRNADisease
load("circ2Disease.dps.RData")##	circRNADisease
load("circFB_dps.RData")
load("circad_use.RData")
data_use = unique(data.frame(circRNA = c(circFB_dps$circRNA, circ2Disease.dps$circRNA.id, circ2Disease.final$ncRNA.id, circad_use$circRNA_ID), Pathway = c(circFB_dps$Pathway, circ2Disease.dps$Pathway, circ2Disease.final$Pathway, circad_use$Pathway)))
length(which(data_use$circRNA %in% circRNA))
save(data_use, file = "circRNA2DiseasePathway.RData")


##	CoExpression
expr_circ = read.table("Healthy_circRNAs.txt", sep = "\t", header = TRUE)
exoRBase_circ = read.table("circRNAs_anno.txt", sep = "\t", header = TRUE)
rns = exoRBase_circ[match(expr_circ$circID, exoRBase_circ$circID), "circBase.ID"]
expr_circ = expr_circ[-which(is.na(rns)), -1]
rownames(expr_circ) = na.omit(rns)
expr_circ = expr_circ[which(rowSums(expr_circ) != 0), ]
expr_m_lnc = read.table("Healthy_longRNAs.txt", sep = "\t", header = TRUE, row.names = 1)
exoRBase_m_lnc = read.table("longRNAs_anno.txt", sep = "\t", header = TRUE)
expr_m = expr_m_lnc[which(exoRBase_m_lnc$Gene.type == "protein coding gene"), ]
rns = NCBI_ID_sub[match(rownames(expr_m), NCBI_ID_sub$Symbol), "GeneID"]
expr_m = expr_m[-which(is.na(rns)), ]##	6915  118
rownames(expr_m) = na.omit(rns)
expr_m = expr_m[which(rowSums(expr_m) != 0), ]##	5972  118
pcc = cor(t(expr_circ), t(expr_m))
save(pcc, file = "circRNA_mRNA_healthy_PCC.RData")
pcc.p = as.data.frame(matrix(NA, nrow = nrow(expr_circ), ncol = nrow(expr_m)), row.names = rownames(expr_circ))
colnames(pcc.p) = rownames(expr_m)
print(dim(pcc.p))
for(i in rownames(pcc.p)){
if(match(i, rownames(pcc.p))%%1000 == 0){print(match(i, rownames(pcc.p)))}
circ.sub = as.numeric(expr_circ[i, ])
pcc.p[i, ] = apply(expr_m, 1, function(x){return(cor.test(circ.sub, x)$p.value)})
}
save(pcc.p, file = "circRNA_mRNA_healthy_PCC_P.RData")
load("circRNA_mRNA_healthy_PCC.RData")
load("circRNA_mRNA_healthy_PCC_P.RData")
pcc.threshold = 0.7
circ.m = list()
for(i in rownames(pcc)[6779:nrow(pcc)]){
circ.m[[i]] = colnames(pcc)[pcc.p[i,] < 0.05 & pcc[i, ] > pcc.threshold]
}
save(circ.m, file = paste0("circRNA_CoExpression_mRNA_pcc", pcc.threshold, ".RData"))
load(paste0("circRNA_CoExpression_mRNA_pcc", pcc.threshold, ".RData"))
table(unlist(lapply(circ.m, length)))

rownames(KEGG) = toupper(rownames(KEGG))
N = 19982##	GENCODE
load(paste0("circRNA_CoExpression_mRNA_pcc", pcc.threshold, ".RData"))
circ.m.HGT = list()
circ.m.pathway = list()
for(circ in names(circ.m)){
if(match(circ, names(circ.m)) %% 1000 ==0){print(match(circ, names(circ.m)))}
if(length(circ.m[[circ]]) != 0){
HGT.res = data.frame()
for(j in 1:nrow(KEGG)){
HGT.res[j, "pathway"] = rownames(KEGG)[j]
pathway.gene = na.omit(as.numeric(KEGG[j, ]))
HGT.res[j, "pathway.gene"] = length(pathway.gene)
x1 = length(intersect(circ.m[[circ]], pathway.gene))
HGT.res[j, "overlap.gene"] = x1
x2 = length(circ.m[[circ]])
x3 = length(pathway.gene)
HGT.res[j, "p.value"] = phyper(x1-1,x2, N-x2, x3, lower.tail=F)
}
HGT.res.order = HGT.res[order(HGT.res$p.value), ]
HGT.res.order$fdr = p.adjust(HGT.res.order$p.value, method="fdr")
circ.m.HGT[[circ]] = HGT.res.order
circ.m.pathway[[circ]] = HGT.res.order[HGT.res.order$p.value < 0.05, "pathway"]
}
}
save(circ.m.HGT, file = paste0("circRNA_CoExpression_mRNA_pcc", pcc.threshold, "_HGT.RData"))
save(circ.m.pathway, file = paste0("circRNA_CoExpression_mRNA_pcc", pcc.threshold, "_pathwayList.RData"))
print(table(unlist(lapply(circ.m.pathway, length))))

load(paste0("circRNA_CoExpression_mRNA_pcc", pcc.threshold, "_pathwayList.RData"))
circ.m.pathway.df = as.data.frame(matrix(NA, nrow = length(pathway), ncol = length(circ.m.pathway)), row.names = pathway)
colnames(circ.m.pathway.df) = names(circ.m.pathway)
for(circ in names(circ.m.pathway)){
circ.m.pathway.df[circ.m.pathway[[circ]], circ] = circ
}
circ.m.pathway.list = list()
for(i in pathway){
circ = as.character(circ.m.pathway.df[i, -which(is.na(circ.m.pathway.df[i,]))])
if(length(circ) != 0){
circ.m.pathway.list[[i]] = circ
}
}
save(circ.m.pathway.list, file = paste0("circRNA_CoExpression_mRNA_pcc", pcc.threshold, "_mRNAList.RData"))
print(table(unlist(lapply(circ.m.pathway.list, length))))

exoRBase_circ = read.table("/XXX/exoRBaseV2/circRNAs_anno.txt", sep = "\t", header = TRUE)
exoRBase_m_lnc = read.table("/XXX/exoRBaseV2/longRNAs_anno.txt", sep = "\t", header = TRUE)
projects = c("BRCA", "CRC", "GBM", "GC", "HCC", "KIRC", "PAAD", "SCLC", "MEL")
TP = toupper(c("Breast cancer", "Colorectal cancer", "GLIOMA", "Gastric cancer", "Hepatocellular carcinoma", "Renal cell carcinoma", "Pancreatic cancer", "SMALL CELL LUNG CANCER", "MELANOMA"))
names(TP) = projects
for(project in projects){
expr_circ = read.table(paste0(filepath, project, "_circRNAs.txt"), header = TRUE, sep = "\t", quote = "", fill = TRUE)
rns = exoRBase_circ[match(expr_circ$circID, exoRBase_circ$circID), "circBase.ID"]
expr_circ = expr_circ[-which(is.na(rns)), -1]
expr_circ = as.data.frame(apply(expr_circ, 2, as.numeric))
expr_circ = aggregate(expr_circ, by = list(na.omit(rns)), FUN = mean)
rownames(expr_circ) = expr_circ[,1]
expr_circ = expr_circ[, -1]
expr_circ = expr_circ[which(rowSums(expr_circ) != 0), ]
expr_circ = log2(expr_circ+1)
expr_m_lnc = read.table(paste0(filepath, project, "_longRNAs.txt"), header = TRUE, sep = "\t", quote = "", fill = TRUE)
expr_m = expr_m_lnc[which(exoRBase_m_lnc$Gene.type == "protein coding gene"), ]
rns = NCBI_ID_sub[match(expr_m[,1], NCBI_ID_sub$Symbol), "GeneID"]
expr_m = expr_m[-which(is.na(rns)), -1]
expr_m = as.data.frame(apply(expr_m, 2, as.numeric))
expr_m = aggregate(expr_m, by = list(na.omit(rns)), FUN = mean)
rownames(expr_m) = expr_m[,1]
expr_m = expr_m[, -1]
expr_m = expr_m[which(rowSums(expr_m) != 0), ]
expr_m = log2(expr_m+1)
pcc = cor(t(expr_circ), t(expr_m))
save(pcc, file = paste0(project, "_circRNA_mRNA_PCC.RData"))
pcc.p = as.data.frame(matrix(NA, nrow = nrow(expr_circ), ncol = nrow(expr_m)), row.names = rownames(expr_circ))
colnames(pcc.p) = rownames(expr_m)
print(dim(pcc.p))
for(i in rownames(pcc.p)){
if(match(i, rownames(pcc.p))%%1000 == 0){print(match(i, rownames(pcc.p)))}
circ.sub = as.numeric(expr_circ[i, ])
pcc.p[i, ] = apply(expr_m, 1, function(x){return(cor.test(circ.sub, x)$p.value)})
}
save(pcc.p, file = paste0(project, "_circRNA_mRNA_PCC_P.RData"))
}

for(project in projects[-1]){
print(project)
load(file = paste0(project, "_circRNA_mRNA_PCC.RData"))
load(file = paste0(project, "_circRNA_mRNA_PCC_P.RData"))
for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
circ.m = list()
for(i in rownames(pcc)){
circ.m[[i]] = colnames(pcc)[pcc.p[i,] < 0.05 & pcc[i, ] > pcc.threshold]
}
save(circ.m, file = paste0(project, "_circRNA_mRNA_PCC_", pcc.threshold, ".RData"))
}
}

table(unlist(lapply(circ.m, length)))
ps aux | grep run_pcc_filter.R


rownames(KEGG) = toupper(rownames(KEGG))
N = 19982##	GENCODE
for(project in projects[-1]){
print(project)
for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
load(paste0(project, "_circRNA_mRNA_PCC_", pcc.threshold, ".RData"))
circ.m.HGT = list()
circ.m.pathway = list()
for(circ in names(circ.m)){
if(length(circ.m[[circ]]) != 0){
HGT.res = data.frame()
for(j in 1:nrow(KEGG)){
HGT.res[j, "pathway"] = rownames(KEGG)[j]
pathway.gene = na.omit(as.numeric(KEGG[j, ]))
HGT.res[j, "pathway.gene"] = length(pathway.gene)
x1 = length(intersect(circ.m[[circ]], pathway.gene))
HGT.res[j, "overlap.gene"] = x1
x2 = length(circ.m[[circ]])
x3 = length(pathway.gene)
HGT.res[j, "p.value"] = phyper(x1-1,x2, N-x2, x3, lower.tail=F)
}
HGT.res.order = HGT.res[order(HGT.res$p.value), ]
HGT.res.order$fdr = p.adjust(HGT.res.order$p.value, method="fdr")
circ.m.HGT[[circ]] = HGT.res.order
circ.m.pathway[[circ]] = HGT.res.order[HGT.res.order$p.value < 0.05, "pathway"]
}
}
save(circ.m.HGT, file = paste0(project, "_circRNA_mRNA_PCC_", pcc.threshold, "_HGT.RData"))
save(circ.m.pathway, file = paste0(project, "_circRNA_mRNA_PCC_", pcc.threshold, "_pathwayList.RData"))
}
}
ps aux | grep run_HGT.R

for(project in projects){
print(project)
for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
circ.m.pathway = list()
load(paste0(project, "_circRNA_mRNA_PCC_", pcc.threshold, "_HGT.RData"))
for(circ in names(circ.m.HGT)){
HGT.res.order = circ.m.HGT[[circ]]
circ.m.pathway[[circ]] = HGT.res.order[HGT.res.order$p.value < 0.05, "pathway"]
}
save(circ.m.pathway, file = paste0(project, "_circRNA_mRNA_PCC_", pcc.threshold, "_pathwayList.RData"))
}
}

for(project in projects){
print(project)
for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
load(paste0(project, "_circRNA_mRNA_PCC_", pcc.threshold, "_pathwayList.RData"))
circ.neighbor.m.pathway.df = as.data.frame(matrix(NA, nrow = length(pathway), ncol = length(circ.m.pathway)), row.names = pathway)
colnames(circ.neighbor.m.pathway.df) = names(circ.m.pathway)
for(circ in names(circ.m.pathway)){
circ.neighbor.m.pathway.df[circ.m.pathway[[circ]], circ] = circ
}
circ.neighbor.m.pathway.list = list()
for(i in pathway){
circs = as.character(circ.neighbor.m.pathway.df[i, -which(is.na(circ.neighbor.m.pathway.df[i,]))])
if(length(circs) != 0){
circ.neighbor.m.pathway.list[[i]] = circs
}
}
circ.m.pathway.list = circ.neighbor.m.pathway.list
save(circ.m.pathway.list, file = paste0(project, "_circRNA_mRNA_PCC_", pcc.threshold, "_mRNAList.RData"))
print(table(unlist(lapply(circ.m.pathway.list, length))))
}
}


##starbase	circRNA-RBP
stb_RBP_circ = read.table("RBP-circRNA-all.txt", header = TRUE, sep = "\t")
stb_RBP_circ = unique(stb_RBP_circ[,c("RBP","geneID","geneName")])
length(which(unique(stb_RBP_circ$RBP) %in% NCBI_ID_sub$Symbol))
unique(stb_RBP_circ$RBP)[-which(unique(stb_RBP_circ$RBP) %in% NCBI_ID_sub$Symbol)]
stb_RBP_circ[which(stb_RBP_circ$RBP == "LIN28"), "RBP"] = "LIN28A"
stb_RBP_circ[which(stb_RBP_circ$RBP == "PAPD5"), "RBP"] = "TENT4B"
stb_RBP_circ[which(stb_RBP_circ$RBP == "RNF219"), "RBP"] = "OBI1"
stb_RBP_circ[which(stb_RBP_circ$RBP == "TROVE2"), "RBP"] = "RO60"
length(which(unique(stb_RBP_circ$geneID) %in% circRNA_info$best_transcript))
length(which(unique(stb_RBP_circ$geneName) %in% circRNA_info$gene_symbol))
stb_RBP_circ$RBP.Entrez = NCBI_ID_sub[match(stb_RBP_circ$RBP, NCBI_ID_sub$Symbol), "GeneID"]
stb_RBP_circ$circBaseID = NCBI_ID_sub[match(stb_RBP_circ$RBP, NCBI_ID_sub$Symbol), "GeneID"]
colnames(stb_RBP_circ)[2:3] = colnames(circRNA_info)[2:3]
temp1 = merge(stb_RBP_circ, circRNA_info[,1:2], by.x = "best_transcript")
temp2 = merge(stb_RBP_circ, circRNA_info[,c(1,3)], by.x = "gene_symbol")
stb_RBP_circ = unique(rbind(temp1[, c("RBP.Entrez", "circRNA_ID")], temp2[, c("RBP.Entrez", "circRNA_ID")]))
save(stb_RBP_circ, file = "stb_RBP_circ.RData")
load("stb_RBP_circ.RData")
CSCD_RBP_circ = read.table("/XXX/CSCD/hg38_common_circrna_rbp.txt", header = FALSE, sep = "\t", quote = "")
CSCD_RBP_circ = unique(CSCD_RBP_circ[, c(1,5)])
CSCD_RBP_circ[,1] = gsub("[|]", "-",CSCD_RBP_circ[,1])
CSCD_RBP_circ$circRNA_ID = circBank[match(CSCD_RBP_circ[,1], circBank$position), "circbaseID"]
circRNA_ID$position = paste(circRNA_ID$chrom, paste(circRNA_ID$start, circRNA_ID$end, sep = "-"), sep = ":")

CSCD_RBP_circ_use = unique(na.omit(CSCD_RBP_circ))
CSCD_RBP_circ_use[,2] = toupper(CSCD_RBP_circ_use[,2])
length(which(unique(CSCD_RBP_circ_use[,2]) %in% NCBI_ID_sub$Symbol))
unique(CSCD_RBP_circ_use[,2])[-which(unique(CSCD_RBP_circ_use[,2]) %in% NCBI_ID_sub$Symbol)]
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "U2AF65"), 2] = "U2AF2"
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "EIF4AIII"), 2] = "EIF4A3"
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "PTB"), 2] = NA
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "FMRP"), 2] = "FMR1"
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "TDP-43"), 2] = "TARDBP"
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "SFRS1"), 2] = "SRSF1"
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "FUS-MUTANT"), 2] = "FUS"
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "LIN28"), 2] = "LIN28A"
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "HUR"), 2] = "ELAVL1"
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "C17ORF85"), 2] = "NCBP3"
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "TNRC6"), 2] = "TNRC6A"
CSCD_RBP_circ_use[which(CSCD_RBP_circ_use[,2] == "C22ORF28"), 2] = "RTCB"
CSCD_RBP_circ_use$RBP.Entrez = NCBI_ID_sub[match(CSCD_RBP_circ_use[,2], NCBI_ID_sub$Symbol), "GeneID"]
CSCD_RBP_circ_use = unique(na.omit(CSCD_RBP_circ_use[, 3:4]))
length(which(unique(CSCD_RBP_circ_use$circRNA_ID) %in% circRNA))##	
save(CSCD_RBP_circ_use, file = "CSCD_RBP_circ_use.RData")

temp1 = stb_RBP_circ[stb_RBP_circ$circRNA_ID  %in% circRNA, ]
temp2 = CSCD_RBP_circ_use[CSCD_RBP_circ_use$circRNA_ID  %in% circRNA, 2:1]
RBP_circ_all = unique(rbind(temp1, temp2))
save(RBP_circ_all, file = "RBP_circ_all.RData")
load("RBP_circ_all.RData")

rownames(KEGG) = toupper(rownames(KEGG))
N = 19982##	GENCODE
circ.m.HGT = list()
circ.m.pathway = list()
for(circ in unique(RBP_circ_all$circRNA_ID)){
if(match(circ, unique(RBP_circ_all$circRNA_ID)) %% 1000 ==0){print(match(circ, unique(RBP_circ_all$circRNA_ID)))}
circ_m = RBP_circ_all[RBP_circ_all$circRNA_ID == circ, "RBP.Entrez"]
if(length(circ_m) != 0){
HGT.res = data.frame()
for(j in 1:nrow(KEGG)){
HGT.res[j, "pathway"] = rownames(KEGG)[j]
pathway.gene = na.omit(as.numeric(KEGG[j, ]))
HGT.res[j, "pathway.gene"] = length(pathway.gene)
x1 = length(intersect(circ_m, pathway.gene))
HGT.res[j, "overlap.gene"] = x1
x2 = length(circ_m)
x3 = length(pathway.gene)
HGT.res[j, "p.value"] = phyper(x1-1,x2, N-x2, x3, lower.tail=F)
}
HGT.res.order = HGT.res[order(HGT.res$p.value), ]
HGT.res.order$fdr = p.adjust(HGT.res.order$p.value, method="fdr")
circ.m.HGT[[circ]] = HGT.res.order
circ.m.pathway[[circ]] = HGT.res.order[HGT.res.order$p.value < 0.05, "pathway"]
}
}
save(circ.m.HGT, file = "circRNA_RBP_mRNA_HGT.RData")
save(circ.m.pathway, file = "circRNA_RBP_mRNA_pathwayList.RData")

load("circRNA_RBP_mRNA_pathwayList.RData")
circ.m.pathway.df = as.data.frame(matrix(NA, nrow = length(pathway), ncol = length(circ.m.pathway)), row.names = pathway)
colnames(circ.m.pathway.df) = names(circ.m.pathway)
for(circ in names(circ.m.pathway)){
circ.m.pathway.df[circ.m.pathway[[circ]], circ] = circ
}
circ.m.pathway.list = list()
for(i in pathway){
circ = as.character(circ.m.pathway.df[i, -which(is.na(circ.m.pathway.df[i,]))])
if(length(circ) != 0){
circ.m.pathway.list[[i]] = circ
}
}
save(circ.m.pathway.list, file = "circRNA_RBP_mRNA_mRNAList.RData")
print(table(unlist(lapply(circ.m.pathway.list, length))))


# database 
setwd(path)
load("circRNA2DiseasePathway.RData")##	data_use

# expression profile
load("exoRBase_healthy_disease_cpm_limma_DEcirc_data_use.RData")

N = length(unique(circRNA_info$circRNA_ID))##	92375
circ_dps.stat = data.frame(disease.pathway = unique(data_use$Pathway), row.names = unique(data_use$Pathway))
for(i in unique(data_use$Pathway)){
circRNA.disease = unique(data_use[data_use$Pathway == i, "circRNA"])
circ_dps.stat[i, "circRNA.disease"] = length(circRNA.disease)
circRNA.pathway = pathway_circRNA_list[[i]]
circ_dps.stat[i, "circRNA.pathway"] = length(circRNA.pathway)
circ_dps.stat[i, "circRNA.overlap"] = length(intersect(circRNA.disease, circRNA.pathway))
circ_dps.stat[i, "pvalue"] = phyper(circ_dps.stat[i, "circRNA.overlap"]-1,circ_dps.stat[i, "circRNA.disease"], N-circ_dps.stat[i, "circRNA.disease"], circ_dps.stat[i, "circRNA.pathway"], lower.tail=F)
}
circ_dps.stat$FDR = p.adjust(circ_dps.stat[, "pvalue"], method = "fdr")
circ_dps.stat = circ_dps.stat[order(-circ_dps.stat$circRNA.disease), ]
write.xlsx(circ_dps.stat, "exoRBase_healthy_disease_cpm_limma_DEcirc_overlap.xlsx")
res.p = data.frame(ncFN = circ_dps.stat$pvalue, row.names = rownames(circ_dps.stat))

for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
load(paste0("/XXX/exoRBaseV2/circRNA_CoExpression_mRNA_pcc", pcc.threshold, "_mRNAList.RData"))
table(unlist(lapply(circ.m.pathway.list, length)))
circ_dps.stat = data.frame(disease.pathway = unique(data_use$Pathway), row.names = unique(data_use$Pathway))
for(i in unique(data_use$Pathway)){
circRNA.disease = unique(data_use[data_use$Pathway == i, "circRNA"])
circ_dps.stat[i, "circRNA.disease"] = length(circRNA.disease)
circRNA.pathway = circ.m.pathway.list[[i]]
circ_dps.stat[i, "circRNA.pathway"] = length(circRNA.pathway)
circ_dps.stat[i, "circRNA.overlap"] = length(intersect(circRNA.disease, circRNA.pathway))
circ_dps.stat[i, "pvalue"] = phyper(circ_dps.stat[i, "circRNA.overlap"]-1,circ_dps.stat[i, "circRNA.disease"], N-circ_dps.stat[i, "circRNA.disease"], circ_dps.stat[i, "circRNA.pathway"], lower.tail=F)
}
circ_dps.stat$FDR = p.adjust(circ_dps.stat[, "pvalue"], method = "fdr")
circ_dps.stat = circ_dps.stat[order(-circ_dps.stat$circRNA.disease), ]
res.p[, paste0("CoExpression_pcc", pcc.threshold)] = circ_dps.stat$pvalue
print(circ_dps.stat[7,])
}

load("circRNA_RBP_mRNA_mRNAList.RData")
table(unlist(lapply(circ.m.pathway.list, length)))
circ_dps.stat = data.frame(disease.pathway = unique(data_use$Pathway), row.names = unique(data_use$Pathway))
for(i in unique(data_use$Pathway)){
circRNA.disease = unique(data_use[data_use$Pathway == i, "circRNA"])
circ_dps.stat[i, "circRNA.disease"] = length(circRNA.disease)
circRNA.pathway = circ.m.pathway.list[[i]]
circ_dps.stat[i, "circRNA.pathway"] = length(circRNA.pathway)
circ_dps.stat[i, "circRNA.overlap"] = length(intersect(circRNA.disease, circRNA.pathway))
circ_dps.stat[i, "pvalue"] = phyper(circ_dps.stat[i, "circRNA.overlap"]-1,circ_dps.stat[i, "circRNA.disease"], N-circ_dps.stat[i, "circRNA.disease"], circ_dps.stat[i, "circRNA.pathway"], lower.tail=F)
}
circ_dps.stat$FDR = p.adjust(circ_dps.stat[, "pvalue"], method = "fdr")
circ_dps.stat = circ_dps.stat[order(-circ_dps.stat$circRNA.disease), ]
res.p$`circRNA-RBP` = circ_dps.stat$pvalue

save(res.p, file = "circRNA2Disease_all_overlap_p.RData")
res.p = res.p[circ_dps.stat[circ_dps.stat$circRNA.disease>20, "disease.pathway"], ]


##	rank
res = data_use
for(i in unique(res$circRNA)){
input = circRNA_GSEA_res[[i]]
res[res$circRNA == i, "ncFN"] = match(res[res$circRNA == i, "Pathway"], input$NAME)
}
for(pcc.threshold in c(0.3, 0.5, 0.7)){
print(pcc.threshold)
load(paste0("/XXX/exoRBaseV2/circRNA_CoExpression_mRNA_pcc", pcc.threshold, "_HGT.RData"))
for(i in intersect(names(circ.m.HGT), unique(res$circRNA))){
input = circ.m.HGT[[i]]
res[res$circRNA == i, paste0("CoExpression_pcc", pcc.threshold)] = match(res[res$circRNA == i, "Pathway"], input$pathway)
}
rm(circ.m.HGT)
}
load("circRNA_RBP_mRNA_HGT.RData")
for(i in intersect(names(circ.m.HGT), unique(res$circRNA))){
input = circ.m.HGT[[i]]
res[res$circRNA == i, "circRNA-RBP"] = match(res[res$circRNA == i, "Pathway"], input$pathway)
}

save(res, file = "DEcircRNA/circRNA2Disease_all_rank.RData")
save(res, file = "circRNA2Disease_all_rank.RData")

# expression profile
load("exoRBase_healthy_disease_cpm_limma_DEcirc.RData")
projects = c("BRCA", "CRC", "GBM", "GC", "HCC", "KIRC", "PAAD", "SCLC")
TP = toupper(c("Breast cancer", "Colorectal cancer", "GLIOMA", "Gastric cancer", "Hepatocellular carcinoma", "Renal cell carcinoma", "Pancreatic cancer", "SMALL CELL LUNG CANCER"))
names(TP) = projects
data_use = data.frame()
for(i in names(res)){
input = subset(res[[i]], abs(logFC)>1&adj.P.Val<0.05)
print(dim(input))
if(nrow(input) != 0){
data_use = rbind(data_use, data.frame(circRNA = rownames(input), Pathway = TP[i]))
}}##	DEcircRNA 
save(data_use, file = "exoRBase_healthy_disease_cpm_limma_DEcirc_data_use.RData")
