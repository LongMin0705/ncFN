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
# miRNA
mRNA = read.table("mRNA.txt")[,1]##	17355
miRNA = read.table("miRNA.txt")[,1]##	1095
miRNA_name = read.xlsx("miRBase_miRNA.xlsx")

load("KEGG.RData")
pathway = toupper(rownames(KEGG))

r=0.2
path = paste0("/XXX/RWR",r,"/res")
setwd(path)

p.thres = data.frame(statistic = c(rep("NOM p-val", 3),rep("FDR q-val", 3)), thres = c(0.001, 0.01, 0.05, 0.05, 0.1, 0.25))
p = 2
load(paste0("miRNA_pathway_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
table(unlist(lapply(miRNA_pathway_list, length)))
load(paste0("pathway_miRNA_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
table(unlist(lapply(pathway_miRNA_list, length)))
load("miRNA_GSEA_res.RData")

path = paste0("/XXX/RWR",r,"/res/",p.thres[p, "statistic"], "_", p.thres[p, "thres"])
# dir.create(path)
setwd(path)
load(paste0(filePath, "miRNA2Disease_all.RData"))##	data_use
load(paste0(filePath, "miRNA2Disease_all_overlap_p.RData"))
# load(paste0(filePath, "miRNA2Disease_all_overlap_fdr.RData"))
res = data.frame(disease.pathway = unique(data_use$Pathway), row.names = unique(data_use$Pathway))
N = length(unique(miRNA_name$ACC))
for(i in res$disease.pathway){
miRNA.disease = unique(data_use[data_use$Pathway == i, "ACC"])
res[i, "miRNA.disease"] = length(miRNA.disease)
miRNA.pathway = pathway_miRNA_list[[i]]
res[i, "miRNA.pathway"] = length(miRNA.pathway)
res[i, "miRNA.overlap"] = length(intersect(miRNA.disease, miRNA.pathway))
res[i, "HGT.pvalue"] = phyper(res[i, "miRNA.overlap"]-1,res[i, "miRNA.disease"], N-res[i, "miRNA.disease"], res[i, "miRNA.pathway"], lower.tail=F)
}
res[res[, "HGT.pvalue"] > 0.05, "disease.pathway"]
res$HGT.fdr = p.adjust(res[, "HGT.pvalue"], method = "fdr")
res[res[, "HGT.fdr"] > 0.05, "disease.pathway"]
res.order = res[order(-res$miRNA.disease), ]
write.xlsx(res.order, "pathway_miRNA_overlap_HGT_db.xlsx")
nrow(res.order)
length(which(res.order$miRNA.disease>=20))
length(which(res.order$miRNA.disease>=20 & res.order$HGT.pvalue<0.05))
res.order = read.xlsx("pathway_miRNA_overlap_HGT_db.xlsx")
# res.order[res.order[, "HGT.fdr"] > 0.05, "disease.pathway"]
res.p$ncFN = res.order[, "HGT.pvalue"]
res.p = res.p[which(res.order$miRNA.disease>=20), ]
res.order = res.order[res.order$miRNA.disease>=20, ]
res.order[res.order[, "HGT.pvalue"] > 0.05, "disease.pathway"]
# res.fdr$ncFN = res.order[, "HGT.fdr"]
compire = list()
for(i in 2:ncol(res.p)){
compire[[i-1]] = c(colnames(res.p)[1], colnames(res.p)[i])
}
data.use = melt(res.p)
data.use$value = -log10(data.use$value)
max(data.use$value)
pdf("pathway_miRNA_overlap_HGT_p_violinPlot_sig_db.pdf")
p = ggplot(data.use, aes(x = variable, y = value))+ geom_violin(aes(fill = variable),scale = "width", alpha=0.7)+geom_boxplot(width = 0.1)+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), text = element_text(size = 12),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 20, face="bold", hjust=0.5))+ scale_y_continuous(name = "-log10(p value)", breaks=seq(0,160,40), labels = seq(0,160,40))+scale_x_discrete(name = "Method")+stat_compare_means(label = "p.signif", method = "t.test", ref.group = "ncFN", label.y = 160)+guides(fill=FALSE)+ scale_fill_brewer(palette = "Set1")
print(p)
dev.off()
load(file = paste0(filePath, "miRNA2Disease_all_rank.RData"))
head(res)
res$ncFN = NA
for(i in intersect(names(miRNA_GSEA_res),unique(res$ACC))){
input = miRNA_GSEA_res[[i]]
input[is.na(input$NES), "NOM p-val"] = 1
input[is.na(input$NES), "NES"] = -1
input[input$NES < 0, "NOM p-val"] = input[input$NES < 0, "NOM p-val"]+1
input.order = input[order(input$`NOM p-val`, -input$NES), ]
res[res$ACC == i, "ncFN"] = match(res[res$ACC == i, "Pathway"], input.order$NAME)
}

for(i in 4:10){
print(wilcox.test(na.omit(res[,3]), na.omit(res[,i]), alternative="less",exact=FALSE,correct=FALSE))
}
#	all signifcant
data.use = melt(res[, 3:ncol(res)])
data.use = na.omit(data.use)
data.use$variable = factor(data.use$variable, levels = colnames(res)[3:ncol(res)])
colnames(data.use)[1] = "method"
pdf("miRNA_method_compare_rank_densityPlot_db.pdf")
ggplot(data.use,aes(value,fill = method, color=method))+geom_density(alpha = 0) + theme_bw()+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), axis.title = element_text(size=16),axis.text=element_text(size=14), plot.title = element_text(size = 18, face="bold", hjust=0.5))+xlab("rank")+ scale_colour_brewer(palette = "Set1")
dev.off()

load(paste0(filePath, "/miRNA2Disease_all_overlap_p.RData"))
temp = read.xlsx(paste0(filePath, "/pathway_miRNA_overlap_HGT.xlsx"))
load(paste0(filePath, "/miRNA2Disease_all_rank.RData"))##	data_use
data_use = res
res = data.frame(disease.pathway = temp$disease.pathway, row.names = temp$disease.pathway)
N = length(unique(miRNA_name$ACC))
for(i in res$disease.pathway){
miRNA.disease = unique(data_use[data_use$Pathway == i, "ACC"])
res[i, "miRNA.disease"] = length(miRNA.disease)
miRNA.pathway = pathway_miRNA_list[[i]]
res[i, "miRNA.pathway"] = length(miRNA.pathway)
res[i, "miRNA.overlap"] = length(intersect(miRNA.disease, miRNA.pathway))
res[i, "HGT.pvalue"] = phyper(res[i, "miRNA.overlap"]-1,res[i, "miRNA.disease"], N-res[i, "miRNA.disease"], res[i, "miRNA.pathway"], lower.tail=F)
}
res[res[, "HGT.pvalue"] > 0.05, "disease.pathway"]
# res$HGT.fdr = p.adjust(res[, "HGT.pvalue"], method = "fdr")
# res[res[, "HGT.fdr"] > 0.05, "disease.pathway"]
res.order = res[order(-res$miRNA.disease), ]
write.xlsx(res.order, "pathway_miRNA_overlap_HGT_expr.xlsx")
nrow(res.order)
length(which(res.order$miRNA.disease>=20))
length(which(res.order$miRNA.disease>=20 & res.order$HGT.pvalue<0.05))
res.order = read.xlsx("pathway_miRNA_overlap_HGT_expr.xlsx")
res.order[res.order[, "HGT.pvalue"] > 0.05, "disease.pathway"]
# res.order[res.order[, "HGT.fdr"] > 0.05, "disease.pathway"]
# res.p$ncFN = res.order[, "HGT.pvalue"]
# res.fdr$ncFN = res.order[, "HGT.fdr"]
compire = list()
for(i in 2:ncol(res.p)){
compire[[i-1]] = c(colnames(res.p)[1], colnames(res.p)[i])
}
data.use = melt(res.p)
data.use$value = -log10(data.use$value)
max(data.use$value)
pdf("pathway_miRNA_overlap_HGT_p_violinPlot_sig_expr.pdf")
ggplot(data.use, aes(x = variable, y = value))+ geom_violin(aes(fill = variable),scale = "width", alpha=0.7)+geom_boxplot(width = 0.1)+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), text = element_text(size = 12),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 20, face="bold", hjust=0.5))+ scale_y_continuous(name = "-log10(p value)", breaks=seq(0,35,10), labels = seq(0,35,10))+scale_x_discrete(name = "Method")+stat_compare_means(label = "p.signif", method = "t.test", ref.group = "TarBase.Weak", label.y = 36)+guides(fill=FALSE)+ scale_fill_brewer(palette = "Set1")
dev.off()
load(paste0(filePath, "/miRNA2Disease_all_rank.RData"))##	res
head(res)
res$ncFN = NA
for(i in intersect(names(miRNA_GSEA_res),unique(res$ACC))){
input = miRNA_GSEA_res[[i]]
input[is.na(input$NES), "NOM p-val"] = 1
input[is.na(input$NES), "NES"] = -1
input[input$NES < 0, "NOM p-val"] = input[input$NES < 0, "NOM p-val"]+1
input.order = input[order(input$`NOM p-val`, -input$NES), ]
res[res$ACC == i, "ncFN"] = match(res[res$ACC == i, "Pathway"], input.order$NAME)
}

for(i in 4:10){
print(wilcox.test(na.omit(res[,3]), na.omit(res[,i]), alternative="less",exact=FALSE,correct=FALSE))
}
#	all signifcant

data.use = melt(res[, 3:ncol(res)])
data.use = na.omit(data.use)
data.use$variable = factor(data.use$variable, levels = colnames(res)[3:ncol(res)])
colnames(data.use)[1] = "method"
pdf("miRNA_method_compare_rank_densityPlot_expr.pdf")
ggplot(data.use,aes(value,fill = method, color=method))+geom_density(alpha = 0) + theme_bw()+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), axis.title = element_text(size=16),axis.text=element_text(size=14), plot.title = element_text(size = 18, face="bold", hjust=0.5))+xlab("rank")+ scale_colour_brewer(palette = "Set1")
dev.off()



# lncRNA
rm(list = ls())
options(stringsAsFactors = FALSE)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(ggsignif)
library(limma)
library(openxlsx)
library(ggpubr)
stat = function(dataset){
  res = data.frame()
  temp = apply(dataset, 2, function(x){
    length(unique(na.omit(x)))
  })
  res = rbind(res,temp)
  temp = apply(dataset, 2, function(x){
    length(na.omit(x))
  })
  res = rbind(res,temp)
  temp = apply(dataset, 2, function(x){
    length(which(is.na(x)))
  })
  res = rbind(res,temp)
  colnames(res) = colnames(dataset)
  rownames(res) = c("unique","All","NA")
  print(res)
}
load("NCBI_ID_sub.RData")
lncRNA = read.table("lncRNA.txt")[,1]
ensemble.lnc = read.table("Homo_sapiens.GRCh38.105.lncRNA.txt")
load("lncipedia_use.RData")
load("KEGG.RData")
pathway = toupper(rownames(KEGG))
p.thres = data.frame(statistic = c(rep("NOM p-val", 3),rep("FDR q-val", 3)), thres = c(0.001, 0.01, 0.05, 0.05, 0.1, 0.25))
p = 2

r = 0.2
path = paste0("/XXX/RWR",r,"/res")
setwd(path)
p.thres[p,]

load(paste0("lncRNA_pathway_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
load(paste0("pathway_lncRNA_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
load("lncRNA_GSEA_res.RData")

path = paste0("/XXX/RWR",r,"/res/",p.thres[p, "statistic"], "_", p.thres[p, "thres"])
setwd(path)

load("lncRNA2Disease_all_overlap_p.RData")##	res.p
res = data.frame(disease.pathway = rownames(res.p), row.names = rownames(res.p))
N = length(unique(ensemble.lnc[,1]))
for(i in res$disease.pathway){
lncRNA.disease = unique(data_use[data_use$Pathway == i, "ncRNA.Ensembl"])
res[i, "lncRNA.disease"] = length(lncRNA.disease)
lncRNA.pathway = pathway_lncRNA_list[[i]]
res[i, "lncRNA.pathway"] = length(lncRNA.pathway)
res[i, "lncRNA.overlap"] = length(intersect(lncRNA.disease, lncRNA.pathway))
res[i, "HGT.pvalue"] = phyper(res[i, "lncRNA.overlap"]-1,res[i, "lncRNA.disease"], N-res[i, "lncRNA.disease"], res[i, "lncRNA.pathway"], lower.tail=F)
}
res[res[, "HGT.pvalue"] > 0.05, "disease.pathway"]
res$HGT.fdr = p.adjust(res[, "HGT.pvalue"], method = "fdr")
res[res[, "HGT.fdr"] > 0.05, "disease.pathway"]
res.order = res[order(-res$lncRNA.disease), ]
write.xlsx(res.order, file = "pathway_lncRNA_overlap_HGT_db.xlsx")
nrow(res.order)
length(which(res.order$lncRNA.disease>=20))
length(which(res.order$lncRNA.disease>=20 & res.order$HGT.pvalue<0.05))
res.order = read.xlsx("pathway_lncRNA_overlap_HGT_db.xlsx")
res.p$ncFN = res.order[, "HGT.pvalue"]
res.p = res.p[,c(1,5:7,2:4)]
compire = list()
for(i in 2:ncol(res.p)){
compire[[i-1]] = c(colnames(res.p)[1], colnames(res.p)[i])
}
data.use = na.omit(melt(res.p))
data.use$value = -log10(data.use$value)
max(data.use$value)
pdf("pathway_lncRNA_overlap_HGT_p_violinPlot_sig_db.pdf")
ggplot(data.use, aes(x = variable, y = value))+ geom_violin(aes(fill = variable),scale = "width", alpha=0.7)+geom_boxplot(width = 0.1)+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), text = element_text(size = 12),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 20, face="bold", hjust=0.5))+ scale_y_continuous(name = "-log10(p value)", breaks=seq(0,70,10), labels = seq(0,70,10))+scale_x_discrete(name = "Method")+stat_compare_means(label = "p.signif", method = "t.test", ref.group = "ncFN", label.y = 67)+guides(fill=FALSE)+ scale_fill_brewer(palette = "Set1")
dev.off()
load(file = paste0(filePath, "lncRNA2Disease_all_rank.RData"))
head(res)
res[res$method=="ncFN", "TPRank"] = NA
for(i in intersect(names(lncRNA_GSEA_res), unique(res$lncRNA))){
input = lncRNA_GSEA_res[[i]]
input[is.na(input$NES), "NOM p-val"] = 1
input[is.na(input$NES), "NES"] = -1
input[input$NES < 0, "NOM p-val"] = input[input$NES < 0, "NOM p-val"]+1
input.order = input[order(input$`NOM p-val`, -input$NES), ]
res[res$lncRNA==i, "TPRank"] = match(res[res$lncRNA == i, "project"], input.order$NAME)
}

for(i in unique(res$method)[-1]){
print(wilcox.test(na.omit(res[res$method=="ncFN",4]), na.omit(res[res$method==i,4]), alternative="less",exact=FALSE,correct=FALSE))
}
#	all signifcant

res$method = factor(res$method, levels = unique(res$method))
pdf("lncRNA_method_compare_rank_densityPlot_db.pdf")
ggplot(res,aes(TPRank,fill = method, color=method))+geom_density(alpha = 0) + theme_bw()+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), axis.title = element_text(size=16),axis.text=element_text(size=14), plot.title = element_text(size = 18, face="bold", hjust=0.5))+ ggtitle("Rank Distribution")+xlab("rank")+ scale_colour_brewer(palette = "Set1")
dev.off()



load(paste0(filePath, "lncRNA2Disease_all_overlap_p.RData"))##	res.p
load(paste0(filePath, "TCGA_FPKM_limma_DElnc_data_use.RData"))##	data_use
head(data_use)
res = data.frame(disease.pathway = rownames(res.p), row.names = rownames(res.p))
N = length(unique(ensemble.lnc[,1]))##	17741
N = 17741
for(i in res$disease.pathway){
lncRNA.disease = unique(data_use[data_use$Pathway == i, "ncRNA.Ensembl"])
res[i, "lncRNA.disease"] = length(lncRNA.disease)
lncRNA.pathway = pathway_lncRNA_list[[i]]
res[i, "lncRNA.pathway"] = length(lncRNA.pathway)
res[i, "lncRNA.overlap"] = length(intersect(lncRNA.disease, lncRNA.pathway))
res[i, "HGT.pvalue"] = phyper(res[i, "lncRNA.overlap"]-1,res[i, "lncRNA.disease"], N-res[i, "lncRNA.disease"], res[i, "lncRNA.pathway"], lower.tail=F)
}
res[res[, "HGT.pvalue"] > 0.05, "disease.pathway"]
res$HGT.fdr = p.adjust(res[, "HGT.pvalue"], method = "fdr")
res[res[, "HGT.fdr"] > 0.05, "disease.pathway"]
res.order = res[order(-res$lncRNA.disease), ]
write.xlsx(res.order, file = "pathway_lncRNA_overlap_HGT_expr.xlsx")
nrow(res.order)
length(which(res.order$lncRNA.disease>=20))
length(which(res.order$lncRNA.disease>=20 & res.order$HGT.pvalue<0.05))
res.order = read.xlsx("pathway_lncRNA_overlap_HGT_expr.xlsx")

res.p$ncFN = res.order[, "HGT.pvalue"]
res.p = res.p[,c(1,5:7,2:4)]

compire = list()
for(i in 2:ncol(res.p)){
compire[[i-1]] = c(colnames(res.p)[1], colnames(res.p)[i])
}
data.use = na.omit(melt(res.p))
data.use$value = -log10(data.use$value)
max(data.use$value)
pdf("pathway_lncRNA_overlap_HGT_p_violinPlot_sig_expr.pdf")
ggplot(data.use, aes(x = variable, y = value))+ geom_violin(aes(fill = variable),scale = "width", alpha=0.7)+geom_boxplot(width = 0.1)+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), text = element_text(size = 12),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 20, face="bold", hjust=0.5))+ scale_y_continuous(name = "-log10(p value)", breaks=seq(0,20,5), labels = seq(0,20,5))+scale_x_discrete(name = "Method")+stat_compare_means(label = "p.signif", method = "t.test", ref.group = "ncFN", label.y = 19)+guides(fill=FALSE)+ scale_fill_brewer(palette = "Set1")
dev.off()
load(file = paste0(filePath, "lncRNA2Disease_all_rank.RData"))
head(res)
res[res$method=="ncFN", "TPRank"] = NA
for(i in intersect(names(lncRNA_GSEA_res), unique(res$lncRNA))){
input = lncRNA_GSEA_res[[i]]
input[is.na(input$NES), "NOM p-val"] = 1
input[is.na(input$NES), "NES"] = -1
input[input$NES < 0, "NOM p-val"] = input[input$NES < 0, "NOM p-val"]+1
input.order = input[order(input$`NOM p-val`, -input$NES), ]
res[res$lncRNA==i, "TPRank"] = match(res[res$lncRNA == i, "project"], input.order$NAME)
}

for(i in unique(res$method)[-1]){
print(wilcox.test(na.omit(res[res$method=="ncFN",4]), na.omit(res[res$method==i,4]), alternative="less",exact=FALSE,correct=FALSE))
}
#	all signifcant

res$method = factor(res$method, levels = unique(res$method))
pdf("lncRNA_method_compare_rank_densityPlot_expr.pdf")
ggplot(res,aes(TPRank,fill = method, color=method))+geom_density(alpha = 0) + theme_bw()+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), axis.title = element_text(size=16),axis.text=element_text(size=14), plot.title = element_text(size = 18, face="bold", hjust=0.5))+ ggtitle("Rank Distribution")+xlab("rank")+ scale_colour_brewer(palette = "Set1")
dev.off()


# circRNA
rm(list = ls())
options(stringsAsFactors = FALSE)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(limma)
library(data.table)
library(RColorBrewer)
library(edgeR)
circBank = read.table("circBank_circrna_annotation.txt", header = TRUE, sep = "\t")
circRNA = read.table("circRNA.txt")[,1]##	7958
circID_to_name <- read.table("hg19_circID_to_name_circBase.txt", header = T, sep = "\t")
colnames(circID_to_name) = c("circRNA_ID", "circRNA_name")
circRNA_ID <- read.table("hsa_hg19_circRNA.txt", header = T, sep = "\t")
circRNA_info <- circRNA_ID[, c("circRNA.ID", "best.transcript", "gene.symbol")]
colnames(circRNA_info) <- c("circRNA_ID", "best_transcript", "gene_symbol")
table(nchar(circRNA_info$circRNA_ID))
load("NCBI_ID_sub.RData")
load("KEGG.RData")
pathway = toupper(rownames(KEGG))


files = dir()
p.thres = data.frame(statistic = c(rep("NOM p-val", 3),rep("FDR q-val", 3)), thres = c(0.001, 0.01, 0.05, 0.05, 0.1, 0.25))
p = 1

print(p)
circRNA_pathway_list = list()
for(i in files){
setwd(paste("/XXX/GSEA_input/GSEA_out", i, sep="/"))
name = unlist(strsplit(i, "[.]"))
gseaRes = as.data.frame(fread(paste0("gsea_report_for_na_pos_", name[3], ".xls")))
length(which(gseaRes[, p.thres[p, "statistic"]] < p.thres[p, "thres"]))
circRNA_pathway_list[[name[1]]] = gseaRes[gseaRes[, p.thres[p, "statistic"]] < p.thres[p, "thres"], "NAME"]
}
length(circRNA_pathway_list)
save(circRNA_pathway_list, file = paste0("/XXX/circRNA_pathway_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
table(unlist(lapply(circRNA_pathway_list, length)))

pathway_circRNA = data.frame(row.names = pathway)
for(i in 1:length(circRNA_pathway_list)){
pathway_circRNA[circRNA_pathway_list[[i]], i] = names(circRNA_pathway_list)[i]
}
pathway_circRNA = t(pathway_circRNA)
pathway_circRNA_list = list()
for(i in 1:ncol(pathway_circRNA)){
pathway_circRNA_list[[colnames(pathway_circRNA)[i]]] = pathway_circRNA[, i][!is.na(pathway_circRNA[, i])]
}
save(pathway_circRNA_list, file = paste0("/XXX/pathway_circRNA_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
table(unlist(lapply(pathway_circRNA_list, length)))

r = 0.2
load(paste0("circRNA_pathway_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))
load(paste0("pathway_circRNA_", p.thres[p, "statistic"], "_", p.thres[p, "thres"], ".RData"))

library(data.table)
library(openxlsx)
files = dir()[grep("^hsa_circ", dir())]
circRNA_GSEA_res = list()
for(i in files){
setwd(paste(path, i, sep="/"))
name = unlist(strsplit(i, "[.]"))
pos = as.data.frame(fread(paste0("gsea_report_for_na_pos_", name[3], ".xls")))
pos.order = pos[order(pos$`NOM p-val`, -pos$NES), ]
neg = as.data.frame(fread(paste0("gsea_report_for_na_neg_", name[3], ".xls")))
if(neg[1,1] != "NAME"){
neg.order = neg[order(neg$`NOM p-val`, -neg$NES), ]
input = rbind(pos.order, neg.order)
}else {input = pos.order}
input[is.na(input$NES), "NOM p-val"] = 1
input[is.na(input$NES), "NES"] = -1
input[input$NES < 0, "NOM p-val"] = input[input$NES < 0, "NOM p-val"]+1
input.order = input[order(input$`NOM p-val`, -input$NES), ]
circRNA_GSEA_res[[name[1]]] = input.order
}
table(unlist(lapply(circRNA_GSEA_res, nrow)))
save(circRNA_GSEA_res, file = "circRNA_GSEA_res.RData")
load("circRNA_GSEA_res.RData")

path = paste0("/XXX/RWR",r,"/res/",p.thres[p, "statistic"], "_", p.thres[p, "thres"])
setwd(path)
getwd()

load("circRNA2Disease_all_overlap_p.RData")##	res.p
load("circRNA2DiseasePathway.RData")##	data_use
res = data.frame(disease.pathway = rownames(res.p), row.names = rownames(res.p))
N = length(unique(circRNA_info$circRNA_ID))##	92375
for(i in res$disease.pathway){
circRNA.disease = unique(data_use[data_use$Pathway == i, "circRNA"])
res[i, "circRNA.disease"] = length(circRNA.disease)
circRNA.pathway = pathway_circRNA_list[[i]]
res[i, "circRNA.pathway"] = length(circRNA.pathway)
res[i, "circRNA.overlap"] = length(intersect(circRNA.disease, circRNA.pathway))
res[i, "HGT.pvalue"] = phyper(res[i, "circRNA.overlap"]-1,res[i, "circRNA.disease"], N-res[i, "circRNA.disease"], res[i, "circRNA.pathway"], lower.tail=F)
}
res[res[, "HGT.pvalue"] > 0.05, "disease.pathway"]
res$HGT.fdr = p.adjust(res[, "HGT.pvalue"], method = "fdr")
res[res[, "HGT.fdr"] > 0.05, "disease.pathway"]
res.order = res[order(-res$circRNA.disease), ]
write.xlsx(res, "pathway_circRNA_overlap_HGT_db.xlsx")
nrow(res.order)
length(which(res.order$circRNA.disease>=20))
length(which(res.order$circRNA.disease>=20 & res.order$HGT.pvalue<0.05))
res = read.xlsx("pathway_circRNA_overlap_HGT.xlsx")
res.p$ncFN = res[, "HGT.pvalue"]
res.p = res.p[res[res$circRNA.disease > 20, "disease.pathway"], ]
compire = list()
for(i in 2:ncol(res.p)){
compire[[i-1]] = c(colnames(res.p)[1], colnames(res.p)[i])
}
data.use = na.omit(melt(res.p))
data.use$value = -log10(data.use$value)
max(data.use$value)
pdf("pathway_circRNA_overlap_HGT_p_violinPlot_sig_db.pdf")
ggplot(data.use, aes(x = variable, y = value))+ geom_violin(aes(fill = variable),scale = "width", alpha=0.7)+geom_boxplot(width = 0.1)+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), text = element_text(size = 12),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 20, face="bold", hjust=0.5))+ scale_y_continuous(name = "-log10(p value)", breaks=seq(0,7,1), labels = seq(0,7,1))+scale_x_discrete(name = "Method")+stat_compare_means(label = "p.signif", method = "t.test", ref.group = "ncFN", label.y = 6.9)+guides(fill=FALSE)+ scale_fill_brewer(palette = "Set1")
dev.off()
load(file = paste0(filePath, "circRNA2Disease_all_rank.RData"))
head(res)
res$ncFN = NA
for(i in intersect(unique(res$circRNA),names(circRNA_GSEA_res))){
input = circRNA_GSEA_res[[i]]
input[is.na(input$NES), "NOM p-val"] = 1
input[is.na(input$NES), "NES"] = -1
input[input$NES < 0, "NOM p-val"] = input[input$NES < 0, "NOM p-val"]+1
input.order = input[order(input$`NOM p-val`, -input$NES), ]
res[res$circRNA == i, "ncFN"] = match(res[res$circRNA == i, "Pathway"], input.order$NAME)
}

for(i in 4:7){
print(wilcox.test(na.omit(res[,3]), na.omit(res[,i]), alternative="less",exact=FALSE,correct=FALSE))
}
#	all signifcant

data.use = melt(res[, 3:ncol(res)])
data.use = na.omit(data.use)
data.use$variable = factor(data.use$variable, levels = colnames(res)[3:ncol(res)])
colnames(data.use) = c("method", "TPRank")
pdf("circRNA_method_compare_rank_densityPlot_db.pdf")
ggplot(data.use,aes(TPRank,fill = method, color=method))+geom_density(alpha = 0) + theme_bw()+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), axis.title = element_text(size=16),axis.text=element_text(size=14), plot.title = element_text(size = 18, face="bold", hjust=0.5))+ ggtitle("Rank Distribution")+xlab("rank")+ scale_colour_brewer(palette = "Set1")
dev.off()


load("circRNA2Disease_all_overlap_p.RData")##	res.p
load("exoRBase_healthy_disease_cpm_limma_DEcirc_data_use.RData")

res = data.frame(disease.pathway = rownames(res.p), row.names = rownames(res.p))
N = length(unique(circRNA_info$circRNA_ID))
for(i in res$disease.pathway){
circRNA.disease = unique(data_use[data_use$Pathway == i, "circRNA"])
res[i, "circRNA.disease"] = length(circRNA.disease)
circRNA.pathway = pathway_circRNA_list[[i]]
res[i, "circRNA.pathway"] = length(circRNA.pathway)
res[i, "circRNA.overlap"] = length(intersect(circRNA.disease, circRNA.pathway))
res[i, "HGT.pvalue"] = phyper(res[i, "circRNA.overlap"]-1,res[i, "circRNA.disease"], N-res[i, "circRNA.disease"], res[i, "circRNA.pathway"], lower.tail=F)
}
res[res[, "HGT.pvalue"] > 0.05, "disease.pathway"]
res$HGT.fdr = p.adjust(res[, "HGT.pvalue"], method = "fdr")
res[res[, "HGT.fdr"] > 0.05, "disease.pathway"]
res.order = res[order(-res$circRNA.disease), ]
write.xlsx(res, "pathway_circRNA_overlap_HGT_expr.xlsx")
nrow(res.order)
length(which(res.order$circRNA.disease>=20))
length(which(res.order$circRNA.disease>=20 & res.order$HGT.pvalue<0.05))
res = read.xlsx("pathway_circRNA_overlap_HGT_expr.xlsx")
res.p$ncFN = res[, "HGT.pvalue"]
compire = list()
for(i in 2:ncol(res.p)){
compire[[i-1]] = c(colnames(res.p)[1], colnames(res.p)[i])
}
data.use = na.omit(melt(res.p))
data.use$value = -log10(data.use$value)
max(data.use$value)
pdf("pathway_circRNA_overlap_HGT_p_violinPlot_sig_expr.pdf")
ggplot(data.use, aes(x = variable, y = value))+ geom_violin(aes(fill = variable),scale = "width", alpha=0.7)+geom_boxplot(width = 0.1)+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), text = element_text(size = 12),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 20, face="bold", hjust=0.5))+ scale_y_continuous(name = "-log10(p value)", breaks=seq(0,120,30), labels = seq(0,120,30))+scale_x_discrete(name = "Method")+stat_compare_means(label = "p.signif", method = "t.test", ref.group = "ncFN", label.y = 127)+guides(fill=FALSE)+ scale_fill_brewer(palette = "Set1")
dev.off()
load(file = paste0(filePath, "circRNA2Disease_all_rank.RData"))
head(res)
res$ncFN = NA
for(i in intersect(unique(res$circRNA),names(circRNA_GSEA_res))){
input = circRNA_GSEA_res[[i]]
input[is.na(input$NES), "NOM p-val"] = 1
input[is.na(input$NES), "NES"] = -1
input[input$NES < 0, "NOM p-val"] = input[input$NES < 0, "NOM p-val"]+1
input.order = input[order(input$`NOM p-val`, -input$NES), ]
res[res$circRNA == i, "ncFN"] = match(res[res$circRNA == i, "Pathway"], input.order$NAME)
}

for(i in 4:7){
print(wilcox.test(na.omit(res[,3]), na.omit(res[,i]), alternative="less",exact=FALSE,correct=FALSE))
}
#	all signifcant

data.use = melt(res[, 3:ncol(res)])
data.use = na.omit(data.use)
data.use$variable = factor(data.use$variable, levels = colnames(res)[3:ncol(res)])
colnames(data.use) = c("method", "TPRank")
pdf("circRNA_method_compare_rank_densityPlot_expr.pdf")
ggplot(data.use,aes(TPRank,fill = method, color=method))+geom_density(alpha = 0) + theme_bw()+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), axis.title = element_text(size=16),axis.text=element_text(size=14), plot.title = element_text(size = 18, face="bold", hjust=0.5))+ ggtitle("Rank Distribution")+xlab("rank")+ scale_colour_brewer(palette = "Set1")
dev.off()
