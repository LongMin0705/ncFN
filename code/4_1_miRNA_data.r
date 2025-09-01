rm(list=ls())
options(stringsAsFactors = FALSE)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(patchwork)
library(openxlsx)
library(ggsignif)

# miRNA
NCBI_ID_sub = read.xlsx("NCBI_ID.xlsx")
mRNA = read.table("mRNA.txt")[,1]
miRNA = read.table("miRNA.txt")[,1]
miRNA_name = read.xlsx("miRBase_miRNA.xlsx")

##	1. HMDD
hmdd = read.xlsx("HMDD.xlsx")
load("HMDD_disease_alisa.RData")
length(disease.alias)
disease.alias.df = data.frame()
for(i in names(disease.alias)){
disease.alias.df = rbind(disease.alias.df, cbind(pathway = i, disease = disease.alias[[i]]))
}
hmdd.pri = unique(hmdd[, 2:3])
hmdd.pri = hmdd.pri[grep("hsa", hmdd.pri$mir), ]
hmdd.pri$disease = toupper(hmdd.pri$disease)
hmdd.pri.sub = merge(hmdd.pri, disease.alias.df, by = "disease")
save(hmdd.pri.sub, file = "hmdd.pri.sub.RData")
load("hmdd.pri.sub.RData")

##	2. miR2Disease
m2D = read.table("miR2Disease.txt", sep = "\t", header = FALSE, quote = "", fill = TRUE)
m2D.pri = unique(m2D[, 1:2])
m2D.pri = m2D.pri[grep("^hsa", m2D.pri[,1]), ]
m2D.pri[,2] = toupper(m2D.pri[,2])
colnames(m2D.pri) = c("miRNA", "Disease")
load("/XXX/miR2Disease/disease.alias.m2D.RData")
disease.alias.m2D.df = data.frame()
for(i in names(disease.alias.m2D)){
disease.alias.m2D.df = rbind(disease.alias.m2D.df, data.frame(Pathway = i, Disease = disease.alias.m2D[[i]]))
}
disease.alias.m2D.df = na.omit(disease.alias.m2D.df)
m2D.pri.sub = merge(m2D.pri, disease.alias.m2D.df, by = "Disease")
length(which(unique(m2D.pri.sub$miRNA) %in% miRNA_name$pri.name))
length(which(unique(m2D.pri.sub$miRNA) %in% miRNA_name$name))
save(m2D.pri.sub, file = "m2D.pri.sub.RData")
load("m2D.pri.sub.RData")

##	3. miRCancer
mc = read.table("miRCancerJune2020.txt", sep = "\t", header = TRUE, quote = "", fill = TRUE)
mc.pri = unique(mc[, 1:2])
mc.pri = mc.pri[grep("hsa", mc.pri$mirId), ]
mc.pri[,2] = toupper(mc.pri[,2])
colnames(mc.pri) = c("miRNA", "Disease")
load("/XXX/miRCancer/disease.alias.mc.RData")
disease.alias.mc.df = data.frame()
for(i in names(disease.alias.mc)){
disease.alias.mc.df = rbind(disease.alias.mc.df, data.frame(Pathway = i, Disease = disease.alias.mc[[i]]))
}
disease.alias.mc.df = na.omit(disease.alias.mc.df)
mc.pri.sub = merge(mc.pri, disease.alias.mc.df, by = "Disease")
length(which(unique(mc.pri.sub$miRNA) %in% miRNA_name$pri.name))
length(which(unique(mc.pri.sub$miRNA) %in% miRNA_name$name))
save(mc.pri.sub, file = "mc.pri.sub.RData")
load("mc.pri.sub.RData")

##	4. TCGA
projects = c("BRCA", "KIRC", "BLCA", "COAD", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "STAD", "THCA", "UCEC")
TP = toupper(c("Breast cancer", "Renal cell carcinoma", "Bladder cancer", "Colorectal cancer", "Hepatocellular carcinoma", "Non-small cell lung cancer", "Non-small cell lung cancer", "Pancreatic cancer", "Prostate cancer", "Gastric cancer", "Thyroid cancer", "Endometrial cancer"))
names(TP) = projects
TCGA.pri.all = data.frame()
for(project in projects){
print(project)
demir = read.xlsx(paste("/XXX", project,"deseq2_res.xlsx", sep = "/"), rowNames = TRUE)
TCGA.pri.all = rbind(TCGA.pri.all, data.frame(Disease = project, pre.name = rownames(demir), Pathway = TP[project]))
}
length(which(unique(TCGA.pri.all$pre.name) %in% miRNA_name$pri.name))
length(which(unique(TCGA.pri.all$pre.name) %in% miRNA_name$name))
save(TCGA.pri.all, file = "TCGA.pri.all.RData")
load("TCGA.pri.all.RData")

projects = c("BRCA", "KIRC", "BLCA", "COAD", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "STAD", "THCA", "UCEC")
TP = toupper(c("Breast cancer", "Renal cell carcinoma", "Bladder cancer", "Colorectal cancer", "Hepatocellular carcinoma", "Non-small cell lung cancer", "Non-small cell lung cancer", "Pancreatic cancer", "Prostate cancer", "Gastric cancer", "Thyroid cancer", "Endometrial cancer"))
names(TP) = projects
TCGA.pri.sub = data.frame()
for(project in projects){
print(project)
demir = read.xlsx(paste("/XXX", project,"deseq2_res.xlsx", sep = "/"), rowNames = TRUE)
demir.sub = subset(demir, padj<0.05&abs(log2FoldChange)>1)
TCGA.pri.sub = rbind(TCGA.pri.sub, data.frame(Disease = project, pre.name = rownames(demir.sub), Pathway = TP[project]))
}
save(TCGA.pri.sub, file = "TCGA.pri.sub.RData")
load("TCGA.pri.sub.RData")

##	5. tarbase
tarbase = read.table("TarBase_v8_download.txt", header = TRUE, sep = "\t")
tarbase.hsa = tarbase[tarbase$species=="Homo sapiens",]
table(tarbase.hsa$direct_indirect)
unique(tarbase.hsa[tarbase.hsa$direct_indirect == "DIRECT", "method"])
unique(tarbase.hsa[tarbase.hsa$direct_indirect == "INDIRECT", "method"])

ids = tarbase.hsa[grep("(hsa)", tarbase.hsa$geneId),"geneId"]
for(i in 1:length(ids)){
ids[i] = substr(ids[i], 1, nchar(ids[i])-5)
}
tarbase.hsa[grep("(hsa)", tarbase.hsa$geneId),"geneId"] = ids
ids = tarbase.hsa[grep("(hsa)", tarbase.hsa$geneName),"geneName"]
for(i in 1:length(ids)){
ids[i] = substr(ids[i], 1, nchar(ids[i])-5)
}
tarbase.hsa[grep("(hsa)", tarbase.hsa$geneName),"geneName"] = ids
length(which(tarbase.hsa$geneName %in% NCBI_ID_sub$Symbol))
tarbase.hsa.sub = unique(tarbase.hsa[, c("geneId","geneName","mirna","method","direct_indirect")])
entrez1 = NCBI_ID_sub[match(tarbase.hsa.sub$geneName, NCBI_ID_sub$Symbol), "GeneID"]
length(na.omit(entrez1))
entrez2 = NCBI_ID_sub[match(tarbase.hsa.sub$geneId, NCBI_ID_sub$Ensembl), "GeneID"]
length(na.omit(entrez2))
entrez = entrez1
entrez[which(is.na(entrez))] = entrez2[which(is.na(entrez))]
length(na.omit(entrez))
tarbase.hsa.sub = unique(cbind(Gene = entrez, tarbase.hsa.sub[, 3:5]))
tarbase.hsa.sub = unique(na.omit(tarbase.hsa.sub))
colnames(tarbase.hsa.sub)[2] = "miRNA"
save(tarbase.hsa.sub, file = "tarbase.hsa.sub.RData")
load("tarbase.hsa.sub.RData")

##	6. mirTarBase
mtb = read.xlsx("hsa_MTI_mirTarBase.xlsx")
table(mtb$`Species.(miRNA)`)
table(mtb$Support.Type)
mtb.all = unique(mtb[, c("miRNA", "Target.Gene.(Entrez.Gene.ID)", "Support.Type")])
colnames(mtb.all)[2] = "Gene"
save(mtb.all, file = "mtb.all.RData")
load("mtb.all.RData")

##	7. StarBase
stb = read.table("StarBase_miRNA_mRNA.txt", sep = "\t", header = TRUE)
stb.use = stb[, c(1,4)]
colnames(stb.use) = c("ACC", "Gene")
save(stb.use, file = "stb.use.RData")
load("stb.use.RData")

##	8. mirDB
load("miRDB_multiMir_network.RData")
miRDB.use = unique(miRDB[, c("mature_mirna_acc","target_entrez")])
colnames(miRDB.use) = c("ACC", "Gene")
save(miRDB.use, file = "miRDB.use.RData")
load("miRDB.use.RData")

##	9. TargetScan
load("TargetScan_multiMir_network.RData")
tgs.use = unique(tgs[, c("mature_mirna_acc","target_entrez")])
colnames(tgs.use) = c("ACC", "Gene")
save(tgs.use, file = "tgs.use.RData")
load("tgs.use.RData")

miRNA.pri = unique(c(hmdd.pri.sub$mir, m2D.pri.sub$miRNA, mc.pri.sub$miRNA, TCGA.pri.all$pre.name, tarbase.hsa.sub$miRNA, mtb.all$miRNA))
length(which(miRNA.pri %in% miRNA_name$pri.name))
length(which(miRNA.pri %in% miRNA_name$name))
miRNA.miss = miRNA.pri[-which(miRNA.pri %in% miRNA_name$pri.name | miRNA.pri %in% miRNA_name$name)]
miRNA.version.NA = c()

miRNA.miss.df = data.frame()
for(i in miRNA.miss){
versions = checkMiRNAVersion(i, verbose = F)
temp = miRNA_NameToAccession(i,version = versions)
colnames(temp) = c("pre.name","pre.ACC")
miRNA.miss.df = rbind(miRNA.miss.df, temp)
}
all(grepl("hsa", unique(miRNA.miss.df$pre.name)))
save(miRNA.miss.df, file = "miRNA.miss.df1.RData")
load("miRNA.miss.df1.RData")
head(miRNA.miss.df[grep("MI0", miRNA.miss.df$pre.ACC), ])
miRNA_name[miRNA_name$pri.ACC == "MI0003672", ]
head(miRNA.miss.df[grep("MIMAT", miRNA.miss.df$pre.ACC), ])
miRNA_name[miRNA_name$ACC == "MI0003672", ]

idss = which(is.na(miRNA.miss.df$pre.ACC) == FALSE)
head(miRNA.miss.df[idss, ])
for(i in idss){
if(!miRNA.miss.df[i, "pre.ACC"] %in% union(miRNA_name$pri.ACC, miRNA_name$ACC)){
miRNA.version.NA = c(miRNA.version.NA, miRNA.miss.df[i, "pre.name"])
miRNA.miss.df[i, "pre.ACC"] = NA
}
}
length(grep("^MIMAT", miRNA.miss.df$pre.ACC))
length(grep("^MI0", miRNA.miss.df$pre.ACC))

idss = intersect(which(is.na(miRNA.miss.df$pre.ACC)), grep("mir", miRNA.miss.df$pre.name))
for(i in idss){
a = miRNA.miss.df[i, "pre.name"]
miRNA.miss.df[i, "pre.ACC"] = miRNA_name[match(gsub("mir", "miR", a), miRNA_name$name), "ACC"]
}
head(na.omit(miRNA.miss.df[idss,]))
tail(na.omit(miRNA.miss.df[idss,]))

length(grep("let", miRNA.miss.df$pre.name))
head(miRNA.miss.df)
idss = intersect(which(is.na(miRNA.miss.df$pre.ACC)), grep("miR", miRNA.miss.df$pre.name))
for(i in idss){
a = miRNA.miss.df[i, "pre.name"]
b = gsub("miR", "mir", a)
if(b %in% miRNA_name$pri.name){
miRNA.miss.df = rbind(miRNA.miss.df, data.frame(pre.name = a, pre.ACC = miRNA_name[miRNA_name$pri.name==b, "ACC"]))
}
}
length(unique(miRNA.miss.df[926:977, "pre.name"]))
miRNA.miss.df = miRNA.miss.df[-which(miRNA.miss.df$pre.name %in% unique(miRNA.miss.df[926:977, "pre.name"]) & is.na(miRNA.miss.df$pre.ACC)), ]

idss = which(is.na(miRNA.miss.df$pre.ACC))
head(miRNA.miss.df[idss,])
for(i in idss){
a = miRNA.miss.df[i, "pre.name"]
if(grepl("miR", a)){
b = gsub("miR", "mir", a)
}else if(grepl("mir", a)){
b = gsub("mir", "miR", a)
}
versions = checkMiRNAVersion(b, verbose = F)
if(!is.na(miRNA_NameToAccession(b,version = versions)[,2])){
miRNA.miss.df = rbind(miRNA.miss.df, data.frame(pre.name = a, pre.ACC = miRNA_NameToAccession(b,version = versions)[,2]))
}
}
head(miRNA.miss.df[951:960, ])
length(unique(miRNA.miss.df[951:1090, "pre.name"]))
save(miRNA.miss.df, file = "miRNA.miss.df2.RData")
load("miRNA.miss.df2.RData")

for(i in 951:1090){
if(!miRNA.miss.df[i, "pre.ACC"] %in% union(miRNA_name$pri.ACC, miRNA_name$ACC)){
miRNA.version.NA = c(miRNA.version.NA, miRNA.miss.df[i, "pre.name"])
miRNA.miss.df[i, "pre.ACC"] = NA
}
}
idss = which(is.na(miRNA.miss.df[, "pre.ACC"]))
miRNA.miss.df = miRNA.miss.df[-idss[idss>950], ]
length(grep("^MIMAT", miRNA.miss.df[951:1028, "pre.ACC"]))
length(grep("^MI0", miRNA.miss.df[951:1028, "pre.ACC"]))
miRNA.miss.df = miRNA.miss.df[-which(miRNA.miss.df$pre.name %in% unique(miRNA.miss.df[951:1028, "pre.name"]) & is.na(miRNA.miss.df$pre.ACC)), ]

load("miRNA.miss.all.id.df.RData")
idss = which(is.na(miRNA.miss.df$pre.ACC))
miRNA.miss.df[idss, "pre.name"]
head(miRNA.miss.df[idss,])
a = unique(miRNA.miss.all.id.df[miRNA.miss.all.id.df$pre.name %in% miRNA.miss.df[idss, "pre.name"], 1:2])
miRNA.miss.df[miRNA.miss.df$pre.name %in% a$pre.name, ]
for(i in 1:nrow(a)){
miRNA.miss.df[miRNA.miss.df$pre.name == a[i, "pre.name"], "pre.ACC"] = a[i, "pri.ACC"]
}
idss = which(is.na(miRNA.miss.df$pre.ACC))
a = unique(miRNA.miss.all.id.df[miRNA.miss.all.id.df$pre.name %in% miRNA.miss.df[idss, "pre.name"], 1:2])
save(miRNA.miss.df, file = "miRNA.miss.df3.RData")
load("miRNA.miss.df3.RData")

miRNA.version.NA = unique(miRNA.version.NA)
length(intersect(miRNA.miss.df$pre.name, miRNA.version.NA))
a = miRNA.miss.df[miRNA.miss.df$pre.name %in% miRNA.version.NA, ]
write.xlsx(a, "miRNA.version.NA.xlsx")

miRNA.miss.mature.df = miRNA.miss.df[grep("MIMAT", miRNA.miss.df$pre.ACC), ]
colnames(miRNA.miss.mature.df) = c("pre.name", "ACC")
miRNA.miss.mature.df = merge(miRNA.miss.mature.df, miRNA_name, by = "ACC")
miRNA.miss.mature.df = miRNA.miss.mature.df[, c(2,3,4,1,5)]

miRNA.miss.pri.df = miRNA.miss.df[grep("MI0", miRNA.miss.df$pre.ACC), ]
colnames(miRNA.miss.pri.df) = c("pre.name", "pri.ACC")
miRNA.miss.pri.df = merge(miRNA.miss.pri.df, miRNA_name, by = "pri.ACC")
miRNA.miss.pri.df = miRNA.miss.pri.df[, c(2,1,3,4,5)]

miRNA.miss.all.df = unique(rbind(miRNA.miss.mature.df, miRNA.miss.pri.df))
save(miRNA.miss.all.df, file = "miRNA.miss.all.df.RData")
load("miRNA.miss.all.df.RData")


pri.name = miRNA_name$pri.name
name = tolower(miRNA_name$name)
all(pri.name != name)
head(miRNA_name[which(pri.name == name), "pri.ACC"])
miRNA_name_t = unique(data.frame(ACC = c(miRNA_name$pri.ACC, miRNA_name$ACC), name = c(miRNA_name$pri.name, tolower(miRNA_name$name))))
length(unique(miRNA_name[which(pri.name == name), "pri.ACC"]))
miRNA_name_t = miRNA_name_t[-match(unique(miRNA_name[which(pri.name == name), "pri.ACC"]), miRNA_name_t$ACC), ]

head(miRNA.pri)
miRNA.pri.t = data.frame(raw = miRNA.pri, lower = tolower(miRNA.pri), upper = gsub("-mir-", "-miR-", miRNA.pri))
all(grepl("^hsa-", miRNA.pri.t$raw))
grep("^hsa-miR-", miRNA.pri.t$lower)
grep("^hsa-mir-", miRNA.pri.t$upper)
length(which(miRNA.pri.t$lower %in% miRNA_name_t$name))
miRNA.pri.t$pre.ACC = miRNA_name_t[match(miRNA.pri.t$lower, miRNA_name_t$name), "ACC"]
head(miRNA.pri.t)
idss = which(is.na(miRNA.pri.t$pre.ACC))##	372

miRNA.pri.t[i, ]
for(i in idss){
versions = checkMiRNAVersion(miRNA.pri.t[i, "raw"], verbose = F)
miRNA.pri.t[i, "pre.ACC"] = miRNA_NameToAccession(miRNA.pri.t[i, "raw"],version = versions)[,2]
if(is.na(miRNA.pri.t[i, "pre.ACC"])){
if(miRNA.pri.t[i, "raw"] != miRNA.pri.t[i, "lower"]){
versions = checkMiRNAVersion(miRNA.pri.t[i, "lower"], verbose = F)
miRNA.pri.t[i, "pre.ACC"] = miRNA_NameToAccession(miRNA.pri.t[i, "lower"],version = versions)[,2]
}else if(miRNA.pri.t[i, "raw"] != miRNA.pri.t[i, "upper"]){
versions = checkMiRNAVersion(miRNA.pri.t[i, "upper"], verbose = F)
miRNA.pri.t[i, "pre.ACC"] = miRNA_NameToAccession(miRNA.pri.t[i, "upper"],version = versions)[,2]
}
}
}
length(which(miRNA.pri.t$pre.ACC %in% miRNA_name_t$ACC))
miRNA.pri.t[which(miRNA.pri.t$pre.ACC %in% miRNA_name_t$ACC == FALSE), "pre.ACC"] = NA


data_t = hmdd.pri.sub
colnames(data_t) = c("Disease", "pre.name", "Pathway")
head(data_t)
data_t$ACC = NA
idss = which(data_t$pre.name %in% miRNA_name$pri.name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
data_ts = data.frame()
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA_name[miRNA_name$pri.name == a, "ACC"]))
}
head(data_ts)
idss = which(data_t$pre.name %in% miRNA_name$name)
length(idss)
idss = which(data_t$pre.name %in% miRNA.miss.all.df$pre.name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA.miss.all.df[miRNA.miss.all.df$pre.name == a, "ACC"]))
}
head(data_ts)
hmdd.pri.sub.t = data_ts
save(hmdd.pri.sub.t, file = "hmdd.pri.sub.t.RData")
load("hmdd.pri.sub.t.RData")


data_t = m2D.pri.sub
colnames(data_t) = c("Disease", "pre.name", "Pathway")
head(data_t)
data_t$ACC = NA
idss = which(data_t$pre.name %in% miRNA_name$pri.name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
data_ts = data.frame()
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA_name[miRNA_name$pri.name == a, "ACC"]))
}
head(data_ts)
idss = which(data_t$pre.name %in% miRNA_name$name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA_name[miRNA_name$name == a, "ACC"]))
}
head(data_ts)
idss = which(data_t$pre.name %in% miRNA.miss.all.df$pre.name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA.miss.all.df[miRNA.miss.all.df$pre.name == a, "ACC"]))
}
head(data_ts)
m2D.pri.sub.t = data_ts
save(m2D.pri.sub.t, file = "m2D.pri.sub.t.RData")
load("m2D.pri.sub.t.RData")

a = m2D.pri.sub.t$pre.name %in% union(miRNA_name$pri.name, miRNA_name$name)
table(a)
b = unique(m2D.pri.sub.t[!a, c(2,4)])
table(b$pre.name %in% union(miRNA_name$pri.name, miRNA_name$name))
temp.df = data.frame()
for(i in unique(b$pre.name)){
versions = checkMiRNAVersion(i, verbose = F)
temp = miRNA_NameToAccession(i,version = versions)
colnames(temp) = c("pre.name","pre.ACC")
temp.df = rbind(temp.df, temp)
}
table(unique(m2D.pri.sub.t$pre.name) %in% union(miRNA_name$pri.name, miRNA_name$name))


data_t = mc.pri.sub
colnames(data_t) = c("Disease", "pre.name", "Pathway")
head(data_t)
data_t$ACC = NA
idss = which(data_t$pre.name %in% miRNA_name$pri.name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
data_ts = data.frame()
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA_name[miRNA_name$pri.name == a, "ACC"]))
}
head(data_ts)
idss = which(data_t$pre.name %in% miRNA_name$name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA_name[miRNA_name$name == a, "ACC"]))
}
head(data_ts)
idss = which(data_t$pre.name %in% miRNA.miss.all.df$pre.name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA.miss.all.df[miRNA.miss.all.df$pre.name == a, "ACC"]))
}
head(data_ts)
mc.pri.sub.t = data_ts
save(mc.pri.sub.t, file = "mc.pri.sub.t.RData")
load("mc.pri.sub.t.RData")

# summary
db = "HMDD_miR2Disease_miRCancer"
filePath = "HMDD_miR2Disease_miRCancer/"
data_use = unique(rbind(rbind(hmdd.pri.sub.t, m2D.pri.sub.t), mc.pri.sub.t))
save(data_use, file = paste0(filePath, "miRNA2Disease_all.RData"))
load(paste0(filePath, "miRNA2Disease_all.RData"))


data_t = TCGA.pri.all
head(data_t)
data_t$ACC = NA
idss = which(data_t$pre.name %in% miRNA_name$pri.name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
data_ts = data.frame()
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA_name[miRNA_name$pri.name == a, "ACC"]))
}
head(data_ts)
idss = which(data_t$pre.name %in% miRNA_name$name)
length(idss)
idss = which(data_t$pre.name %in% miRNA.miss.all.df$pre.name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA.miss.all.df[miRNA.miss.all.df$pre.name == a, "ACC"]))
}
head(data_ts)
TCGA.pri.all.t = data_ts
save(TCGA.pri.all.t, file = "TCGA.pri.all.t.RData")
load("TCGA.pri.all.t.RData")


data_t = TCGA.pri.sub
head(data_t)
data_t$ACC = NA
idss = which(data_t$pre.name %in% miRNA_name$pri.name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
data_ts = data.frame()
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA_name[miRNA_name$pri.name == a, "ACC"]))
}
head(data_ts)
idss = which(data_t$pre.name %in% miRNA_name$name)
length(idss)
idss = which(data_t$pre.name %in% miRNA.miss.all.df$pre.name)
length(idss)
length(unique(data_t[idss, "pre.name"]))
for(i in idss){
a = data_t[i, "pre.name"]
data_ts = rbind(data_ts, cbind(data_t[i, -4], ACC = miRNA.miss.all.df[miRNA.miss.all.df$pre.name == a, "ACC"]))
}
head(data_ts)
TCGA.pri.sub.t = data_ts
save(TCGA.pri.sub.t, file = "TCGA.pri.sub.t.RData")
load("TCGA.pri.sub.t.RData")

load("tarbase.hsa.sub.RData")
length(grepl("^hsa-miR", tarbase.hsa.sub$miRNA) | grepl("^hsa-let", tarbase.hsa.sub$miRNA))	##	all
tarbase.hsa.sub$ACC = miRNA_name[match(tarbase.hsa.sub$miRNA, miRNA_name$name), "ACC"]
idss = which(is.na(tarbase.hsa.sub$ACC))
length(unique(tarbase.hsa.sub[idss, "miRNA"]))
data_use = unique(miRNA.miss.all.df[miRNA.miss.all.df$pre.name %in% unique(tarbase.hsa.sub[idss, "miRNA"]), c(1,4)])
for(i in data_use$pre.name){
tarbase.hsa.sub[tarbase.hsa.sub$miRNA == i, "ACC"] = data_use[data_use$pre.name==i, "ACC"]
}
idss = which(is.na(tarbase.hsa.sub$ACC))
length(unique(tarbase.hsa.sub[idss, "miRNA"]))
tarbase.use = unique(na.omit(tarbase.hsa.sub[, c(5,1,3,4)]))
save(tarbase.use, file = "tarbase.use.RData")
load("tarbase.use.RData")
tarbase.use.all = unique(tarbase.use[, 1:2])
tarbase.use.low = unique(tarbase.use[tarbase.use$method %in% c("Luciferase Reporter Assay", "Western Blot", "qPCR", "Biotin-qPCR", "Northern Blot"), 1:3])
table(tarbase.use.low$method)
tarbase.use.low = unique(tarbase.use.low[, 1:2])
tarbase.use.high = unique(tarbase.use[!tarbase.use$method %in% c("Luciferase Reporter Assay", "Western Blot", "qPCR", "Biotin-qPCR", "Northern Blot"), 1:3])
table(tarbase.use.high$method)
tarbase.use.high = unique(tarbase.use.high[, 1:2])
table(tarbase.use$direct_indirect)
tarbase.use.direct = unique(tarbase.use[tarbase.use$direct_indirect == "DIRECT", 1:2])
tarbase.use.indirect = unique(tarbase.use[tarbase.use$direct_indirect == "INDIRECT", 1:2])


load("mtb.all.RData")
length(grepl("^hsa-miR", mtb.all$miRNA) | grepl("^hsa-let", mtb.all$miRNA))
mtb.all$ACC = miRNA_name[match(mtb.all$miRNA, miRNA_name$name), "ACC"]
idss = which(is.na(mtb.all$ACC))
length(unique(mtb.all[idss, "miRNA"]))
data_use = unique(miRNA.miss.all.df[miRNA.miss.all.df$pre.name %in% unique(mtb.all[idss, "miRNA"]), c(1,4)])
for(i in data_use$pre.name){
mtb.all[mtb.all$miRNA == i, "ACC"] = data_use[data_use$pre.name==i, "ACC"]
}
idss = which(is.na(mtb.all$ACC))
length(unique(mtb.all[idss, "miRNA"]))
mtb.use = unique(na.omit(mtb.all[, c(4,2,3)]))
save(mtb.use, file = "mtb.use.RData")
load("mtb.use.RData")
mtb.use.all = unique(mtb.use[, 1:2])
mtb.use.low = unique(mtb.use[-grep("Weak", mtb.use$Support.Type), 1:2])
mtb.use.high = unique(mtb.use[grep("Weak", mtb.use$Support.Type), 1:2])

##	miRNA-pathway predict
####	RWR.GSEA
load("miRNA_GSEA_res.RData")
KEGG = read.xlsx("kegg_pathway.xlsx", colNames = FALSE, rowNames = TRUE)
rownames(KEGG) = toupper(rownames(KEGG))
pathway = toupper(rownames(KEGG))
N = 19982

data_use = tarbase.use.high
length(unique(data_use$ACC))
MIT.HGT = list()
for(i in unique(data_use$ACC)){
print(i)
target.gene = unique(data_use[data_use$ACC==i, "Gene"])
res = data.frame()
for(j in 1:nrow(KEGG)){
res[j, "pathway"] = rownames(KEGG)[j]
pathway.gene = na.omit(as.numeric(KEGG[j, ]))
res[j, "pathway.gene"] = length(pathway.gene)
x1 = length(intersect(target.gene, pathway.gene))
res[j, "overlap.gene"] = x1
x2 = length(target.gene)
x3 = length(pathway.gene)
res[j, "p.value"] = phyper(x1-1,x2, N-x2, x3, lower.tail=F)
}
res.order = res[order(res$p.value), ]
res.order$fdr = p.adjust(res.order$p.value, method="fdr")
MIT.HGT[[i]] = res.order
}
table(unlist(lapply(MIT.HGT, nrow)))
save(MIT.HGT, file = "tarbase.use.high.HGT.RData")

####	HGT
rm(MIT.HGT)
load("tarbase.use.all.HGT.RData")
load("tarbase.use.low.HGT.RData")
load("tarbase.use.high.HGT.RData")
load("tarbase.use.direct.HGT.RData")
load("tarbase.use.indirect.HGT.RData")
load("mtb.use.all.HGT.RData")
load("mtb.use.low.HGT.RData")
load("mtb.use.high.HGT.RData")
load("stb.use.HGT.RData")
load("miRDB.use.HGT.RData")
load("tgs.use.HGT.RData")
length(MIT.HGT)
pathway.miRNA.df = data.frame(row.names = pathway)
for(i in names(MIT.HGT)){
input = MIT.HGT[[i]]
input.order = input[match(pathway, input$pathway), ]
pathway.miRNA.df[,i] = ifelse(input.order$p.value < 0.05, 1, 0)
}
table(colSums(pathway.miRNA.df))
table(rowSums(pathway.miRNA.df))
pathway.miRNA.list = list()
for(i in rownames(pathway.miRNA.df)){
pathway.miRNA.list[[i]] = colnames(pathway.miRNA.df)[pathway.miRNA.df[i, ]==1]
}
table(unlist(lapply(pathway.miRNA.list, length)))
save(pathway.miRNA.list, file = "tarbase.use.all.pathway.miRNA.list.RData")
save(pathway.miRNA.list, file = "tarbase.use.low.pathway.miRNA.list.RData")
save(pathway.miRNA.list, file = "tarbase.use.high.pathway.miRNA.list.RData")
save(pathway.miRNA.list, file = "tarbase.use.direct.pathway.miRNA.list.RData")
save(pathway.miRNA.list, file = "tarbase.use.indirect.pathway.miRNA.list.RData")
save(pathway.miRNA.list, file = "mtb.use.all.pathway.miRNA.list.RData")
save(pathway.miRNA.list, file = "mtb.use.low.pathway.miRNA.list.RData")
save(pathway.miRNA.list, file = "mtb.use.high.pathway.miRNA.list.RData")
save(pathway.miRNA.list, file = "stb.use.pathway.miRNA.list.RData")
save(pathway.miRNA.list, file = "miRDB.use.pathway.miRNA.list.RData")
save(pathway.miRNA.list, file = "tgs.use.pathway.miRNA.list.RData")


# 1 database
load("miRNA2Disease_all.RData")

# 2 expr
load("TCGA.pri.sub.t.RData")

# p value of intersection
res.p = data.frame(ncFN = rep(0, length(unique(data_use$Pathway))))
dbs = data.frame(method = c("TarBase.Strong", "TarBase.Weak", "miRTarBase.Strong", "miRTarBase.Weak", "StarBase", "miRDB", "TargetScan"), path = c("tarbase.use.low.pathway.miRNA.list.RData", "tarbase.use.high.pathway.miRNA.list.RData", "mtb.use.low.pathway.miRNA.list.RData", "mtb.use.high.pathway.miRNA.list.RData", "stb.use.pathway.miRNA.list.RData", "miRDB.use.pathway.miRNA.list.RData","tgs.use.pathway.miRNA.list.RData"))
N = length(unique(miRNA_name$ACC))
for(n in 1:nrow(dbs)){
load(dbs[n, "path"])
res = data.frame(disease.pathway = unique(data_use$Pathway), row.names = unique(data_use$Pathway))
for(i in res$disease.pathway){
miRNA.disease = unique(data_use[data_use$Pathway == i, "ACC"])
res[i, "miRNA.disease"] = length(miRNA.disease)
miRNA.pathway = pathway.miRNA.list[[i]]
res[i, "miRNA.pathway"] = length(miRNA.pathway)
res[i, "miRNA.overlap"] = length(intersect(miRNA.disease, miRNA.pathway))
res[i, "HGT.pvalue"] = phyper(res[i, "miRNA.overlap"]-1,res[i, "miRNA.disease"], N-res[i, "miRNA.disease"], res[i, "miRNA.pathway"], lower.tail=F)
}
res.p[, dbs[n, "method"]] = res[, "HGT.pvalue"]
rm(pathway.miRNA.list)
}
save(res.p, file = paste0(filePath, "miRNA2Disease_all_overlap_p.RData"))

# target pathway rank
res = unique(data_use[, 3:4])
for(i in names(miRNA_GSEA_res)){
if(i %in% unique(res$ACC)){
input = miRNA_GSEA_res[[i]]
input[is.na(input$NES), "NOM p-val"] = 1
input[is.na(input$NES), "NES"] = -1
input[input$NES < 0, "NOM p-val"] = input[input$NES < 0, "NOM p-val"]+1
input.order = input[order(input$`NOM p-val`, -input$NES), ]
res[res$ACC == i, "ncFN"] = match(res[res$ACC == i, "Pathway"], input.order$NAME)
}
}
dbs = data.frame(method = c("TarBase.Strong", "TarBase.Weak", "miRTarBase.Strong", "miRTarBase.Weak","StarBase", "miRDB", "TargetScan"), path = c("tarbase.use.low.HGT.RData", "tarbase.use.high.HGT.RData", "mtb.use.low.HGT.RData", "mtb.use.high.HGT.RData", "stb.use.HGT.RData", "miRDB.use.HGT.RData","tgs.use.HGT.RData"))
for(n in 1:nrow(dbs)){
rm(MIT.HGT)
load(dbs[n, "path"])
for(i in names(MIT.HGT)){
if(i %in% unique(res$ACC)){
input = MIT.HGT[[i]]
input.order = input[order(input$p.value), ]
res[res$ACC == i, dbs[n, "method"]] = match(res[res$ACC == i, "Pathway"], input.order$pathway)
}
}
}
save(res, file = paste0(filePath, "miRNA2Disease_all_rank.RData"))
