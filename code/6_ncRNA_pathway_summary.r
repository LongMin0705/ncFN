rm(list = ls())
options(stringsAsFactors = FALSE)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(patchwork)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(data.table)
library(openxlsx)

load("KEGG.RData")
pathway = toupper(rownames(KEGG))
lncRNA = read.table("lncRNA.txt", sep = "\t", header = FALSE)[,1]
mRNA = read.table("mRNA.txt", sep = "\t", header = FALSE)[,1]
miRNA = read.table("miRNA.txt", sep = "\t", header = FALSE)[,1]
circRNA = read.table("circRNA.txt", sep = "\t", header = FALSE)[,1]

r = 0.2
p.thres = data.frame(statistic = c(rep("NOM p-val", 3),rep("FDR q-val", 3)), thres = c(0.001, 0.01, 0.05, 0.05, 0.1, 0.25))
for(p in 1:nrow(p.thres)){
print(p.thres[p, ])
path = paste0("/XXX/RWR",r,"/GSEA")
setwd(path)
files = dir()[grep("MIMAT", dir())]
print(length(files))
miRNA_pathway_list = list()
for(i in files){
setwd(paste(path, i, sep="/"))
name = unlist(strsplit(i, "[.]"))
pos = as.data.frame(fread(paste0("gsea_report_for_na_pos_", name[3], ".xls")))
pos.order = pos[order(pos$`NOM p-val`, -pos$NES), ]
neg = as.data.frame(fread(paste0("gsea_report_for_na_neg_", name[3], ".xls")))
miRNA_pathway_list[[name[1]]] = pos[pos[, p.thres[p, "statistic"]] < p.thres[p, "thres"], "NAME"]
}
save(miRNA_pathway_list, file = paste0("/XXX/RWR", r, "/res/", "miRNA_pathway_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
pathway_miRNA = data.frame(row.names = pathway)
for(i in 1:length(miRNA_pathway_list)){
pathway_miRNA[miRNA_pathway_list[[i]], i] = names(miRNA_pathway_list)[i]
}
pathway_miRNA = t(pathway_miRNA)
pathway_miRNA_list = list()
for(i in 1:ncol(pathway_miRNA)){
pathway_miRNA_list[[colnames(pathway_miRNA)[i]]] = pathway_miRNA[, i][!is.na(pathway_miRNA[, i])]
}
save(pathway_miRNA_list, file = paste0("/XXX/RWR", r, "/res/", "pathway_miRNA_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
# lncRNA
setwd(path)
files = dir()[grep("ENSG", dir())]
print(length(files))
lncRNA_pathway_list = list()
for(i in files){
setwd(paste(path, i, sep="/"))
name = unlist(strsplit(i, "[.]"))
pos = as.data.frame(fread(paste0("gsea_report_for_na_pos_", name[3], ".xls")))
pos.order = pos[order(pos$`NOM p-val`, -pos$NES), ]
neg = as.data.frame(fread(paste0("gsea_report_for_na_neg_", name[3], ".xls")))
lncRNA_pathway_list[[name[1]]] = pos[pos[, p.thres[p, "statistic"]] < p.thres[p, "thres"], "NAME"]
}
save(lncRNA_pathway_list, file = paste0("/XXX/RWR", r, "/res/", "lncRNA_pathway_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
pathway_lncRNA = data.frame(row.names = pathway)
for(i in 1:length(lncRNA_pathway_list)){
pathway_lncRNA[lncRNA_pathway_list[[i]], i] = names(lncRNA_pathway_list)[i]
}
pathway_lncRNA = t(pathway_lncRNA)
pathway_lncRNA_list = list()
for(i in 1:ncol(pathway_lncRNA)){
pathway_lncRNA_list[[colnames(pathway_lncRNA)[i]]] = pathway_lncRNA[, i][!is.na(pathway_lncRNA[, i])]
}
save(pathway_lncRNA_list, file = paste0("/XXX/RWR", r, "/res/", "pathway_lncRNA_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
# circRNA
setwd(path)
files = dir()[grep("hsa_circ", dir())]
print(length(files))
circRNA_pathway_list = list()
for(i in files){
setwd(paste(path, i, sep="/"))
name = unlist(strsplit(i, "[.]"))
pos = as.data.frame(fread(paste0("gsea_report_for_na_pos_", name[3], ".xls")))
pos.order = pos[order(pos$`NOM p-val`, -pos$NES), ]
neg = as.data.frame(fread(paste0("gsea_report_for_na_neg_", name[3], ".xls")))
circRNA_pathway_list[[name[1]]] = pos[pos[, p.thres[p, "statistic"]] < p.thres[p, "thres"], "NAME"]
}
save(circRNA_pathway_list, file = paste0("/XXX/RWR", r, "/res/", "circRNA_pathway_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
pathway_circRNA = data.frame(row.names = pathway)
for(i in 1:length(circRNA_pathway_list)){
pathway_circRNA[circRNA_pathway_list[[i]], i] = names(circRNA_pathway_list)[i]
}
pathway_circRNA = t(pathway_circRNA)
pathway_circRNA_list = list()
for(i in 1:ncol(pathway_circRNA)){
pathway_circRNA_list[[colnames(pathway_circRNA)[i]]] = pathway_circRNA[, i][!is.na(pathway_circRNA[, i])]
}
save(pathway_circRNA_list, file = paste0("/XXX/RWR", r, "/res/", "pathway_circRNA_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
}


setwd(path)
load("miRNA_GSEA_res.RData")
miRNA_pathway_list1 = list()
for(i in names(miRNA_GSEA_res)){
miRNA_pathway_list1[[i]] = c(subset(miRNA_GSEA_res[[i]], NES>0&`NOM p-val`<0.01, select=NAME)[,1])
}
quantile(unlist(lapply(miRNA_pathway_list1, length)))

setwd(paste0(path, r, "/res"))
# miRNA
files = dir()[grep("^miRNA_pathway_", dir())]
data_use = list()
for(f in files){
load(f)
thres = substr(f, 15, nchar(f))
thres = substr(thres, 1, nchar(thres)-6)
data_use[[thres]] = as.numeric(lapply(miRNA_pathway_list, length))
}
df <- data.frame(matrix(unlist(data_use), ncol=length(data_use), byrow=FALSE))
names(df) = names(data_use)
stat = data.frame(threshold = names(df),row.names = names(df))
for(i in stat$threshold){
stat[i, "num=0"] = length(which(data_use[[i]]==0))/length(miRNA)
}
stat
data_plot = melt(df)
pdf("miRNA_related_pathway_count_percent_all.pdf",height=10,width=15)
ggplot(data_plot, aes(x = value))+geom_histogram(binwidth=10,color='black',fill='#EF7D55',cex=1,alpha=0.6,boundary=0)+theme_classic(base_size = 18)+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+labs(x="Pathway Count Interval",y='miRNA Count')+theme(axis.text = element_text(color = 'black',size=18),axis.title = element_text(color = 'black',face="bold",size=24),strip.text.x = element_text(face="bold",size=18))+facet_wrap(~variable, nrow=2, scales = "free_y")
dev.off()
# lncRNA
files = dir()[grep("^lncRNA_pathway_", dir())]
data_use = list()
for(f in files){
load(f)
thres = substr(f, 16, nchar(f))
thres = substr(thres, 1, nchar(thres)-6)
data_use[[thres]] = as.numeric(lapply(lncRNA_pathway_list, length))
}
df <- data.frame(matrix(unlist(data_use), ncol=length(data_use), byrow=FALSE))
names(df) = names(data_use)
stat = data.frame(threshold = names(df),row.names = names(df))
for(i in stat$threshold){
stat[i, "num=0"] = length(which(data_use[[i]]==0))/length(lncRNA)
}
stat
data_plot = melt(df)
pdf("lncRNA_related_pathway_count_percent_all.pdf",height=10,width=15)
ggplot(data_plot, aes(x = value))+geom_histogram(binwidth=10,color='black',fill='#6A5D8D',cex=1,alpha=0.6,boundary=0)+theme_classic(base_size = 18)+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+labs(x="Pathway Count Interval",y='lncRNA Count')+theme(axis.text = element_text(color = 'black',size=18),axis.title = element_text(color = 'black',face="bold",size=24),strip.text.x = element_text(face="bold",size=18))+facet_wrap(~variable, nrow=2, scales = "free_y")
dev.off()
# circRNA
files = dir()[grep("^circRNA_pathway_", dir())]
data_use = list()
for(f in files){
load(f)
thres = substr(f, 17, nchar(f))
thres = substr(thres, 1, nchar(thres)-6)
data_use[[thres]] = as.numeric(lapply(circRNA_pathway_list, length))
}
df <- data.frame(matrix(unlist(data_use), ncol=length(data_use), byrow=FALSE))
names(df) = names(data_use)
stat = data.frame(threshold = names(df),row.names = names(df))
for(i in stat$threshold){
stat[i, "num=0"] = length(which(data_use[[i]]==0))/length(circRNA)
}
stat
data_plot = melt(df)
pdf("circRNA_related_pathway_count_percent_all.pdf",height=10,width=15)
ggplot(data_plot, aes(x = value))+geom_histogram(binwidth=10,color='black',fill='#6F976D',cex=1,alpha=0.6,boundary=0)+theme_classic(base_size = 18)+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+labs(x="Pathway Count Interval",y='circRNA Count')+theme(axis.text = element_text(color = 'black',size=18),axis.title = element_text(color = 'black',face="bold",size=24),strip.text.x = element_text(face="bold",size=18))+facet_wrap(~variable, nrow=2, scales = "free_y")
dev.off()


files = dir()[grep("^pathway_miRNA_", dir())]
data_use = list()
for(f in files){
load(f)
thres = substr(f, 15, nchar(f))
thres = substr(thres, 1, nchar(thres)-6)
data_use[[thres]] = as.numeric(lapply(pathway_miRNA_list, length))
}
df <- data.frame(matrix(unlist(data_use), ncol=length(data_use), byrow=FALSE))
names(df) = names(data_use)
stat = data.frame(threshold = names(df),row.names = names(df))
for(i in stat$threshold){
stat[i, "num=0"] = length(which(data_use[[i]]==0))/length(pathway)
}
stat
data_plot = melt(df)
pdf("pathway_related_miRNA_count_percent_all.pdf",height=10,width=15)
ggplot(data_plot, aes(x = value))+geom_histogram(binwidth=50,color='black',fill='#EF7D55',cex=1,alpha=0.6,boundary=0)+theme_classic(base_size = 18)+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+labs(x="miRNA Count Interval",y='Pathway Count')+theme(axis.text = element_text(color = 'black',size=18),axis.title = element_text(color = 'black',face="bold",size=24),strip.text.x = element_text(face="bold",size=18))+facet_wrap(~variable, nrow=2, scales = "free_y")
dev.off()

files = dir()[grep("^pathway_lncRNA_", dir())]
data_use = list()
for(f in files){
load(f)
thres = substr(f, 16, nchar(f))
thres = substr(thres, 1, nchar(thres)-6)
data_use[[thres]] = as.numeric(lapply(pathway_lncRNA_list, length))
}
df <- data.frame(matrix(unlist(data_use), ncol=length(data_use), byrow=FALSE))
names(df) = names(data_use)
stat = data.frame(threshold = names(df),row.names = names(df))
for(i in stat$threshold){
stat[i, "num=0"] = length(which(data_use[[i]]==0))/length(pathway)
}
stat
data_plot = melt(df)
pdf("pathway_related_lncRNA_count_percent_all.pdf",height=10,width=15)
ggplot(data_plot, aes(x = value))+geom_histogram(binwidth=100,color='black',fill='#6A5D8D',cex=1,alpha=0.6,boundary=0)+theme_classic(base_size = 18)+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+labs(x="lncRNA Count Interval",y='Pathway Count')+theme(axis.text = element_text(color = 'black',size=18),axis.title = element_text(color = 'black',face="bold",size=24),strip.text.x = element_text(face="bold",size=18))+facet_wrap(~variable, nrow=2, scales = "free_y")
dev.off()

files = dir()[grep("^pathway_circRNA_", dir())]
data_use = list()
for(f in files){
load(f)
thres = substr(f, 17, nchar(f))
thres = substr(thres, 1, nchar(thres)-6)
data_use[[thres]] = as.numeric(lapply(pathway_circRNA_list, length))
}
df <- data.frame(matrix(unlist(data_use), ncol=length(data_use), byrow=FALSE))
names(df) = names(data_use)
stat = data.frame(threshold = names(df),row.names = names(df))
for(i in stat$threshold){
stat[i, "num=0"] = length(which(data_use[[i]]==0))/length(pathway)
}
stat
data_plot = melt(df)
pdf("pathway_related_circRNA_count_percent_all.pdf",height=10,width=15)
ggplot(data_plot, aes(x = value))+geom_histogram(binwidth=300,color='black',fill='#6F976D',cex=1,alpha=0.6,boundary=0)+theme_classic(base_size = 18)+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+labs(x="circRNA Count Interval",y='Pathway Count')+theme(axis.text = element_text(color = 'black',size=18),axis.title = element_text(color = 'black',face="bold",size=24),strip.text.x = element_text(face="bold",size=18))+facet_wrap(~variable, nrow=2, scales = "free_y")
dev.off()

