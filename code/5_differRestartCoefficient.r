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

lncRNA = read.table("lncRNA.txt", sep = "\t", header = FALSE)[,1]
mRNA = read.table("mRNA.txt", sep = "\t", header = FALSE)[,1]
miRNA = read.table("miRNA.txt", sep = "\t", header = FALSE)[,1]
circRNA = read.table("circRNA.txt", sep = "\t", header = FALSE)[,1]

# miRNA
filePath = "/XXX/"
load(paste0(filePath, "miRNA2Disease_all.RData"))##	data_use
res = unique(data_use[, 3:4])
for(r in seq(0.1,0.9,0.1)){
for(mir in intersect(unique(res$ACC), miRNA)){
load(paste0("/XXX/RWR", r, "/res", "/miRNA_GSEA_res.RData"))
input = miRNA_GSEA_res[[mir]]
input[is.na(input$NES), "NOM p-val"] = 1
input[is.na(input$NES), "NES"] = -1
input[input$NES < 0, "NOM p-val"] = input[input$NES < 0, "NOM p-val"]+1
input.order = input[order(input$`NOM p-val`, -input$NES), ]
res[res$ACC == mir, paste0("R", r)] = match(res[res$ACC == mir, "Pathway"], input.order$NAME)
}
}
save(res, file = "miRNA_differ_R_rank_comparison.RData")
res = na.omit(res)
cor(res[,-c(1,2)])
library(ggcorrplot)
corr <- round(cor(res[,-c(1,2)]), 1)
pdf("miRNA_differ_R_rank_cor.pdf")
ggcorrplot(corr,hc.order = FALSE,type = "upper",outline.color = "white",lab=TRUE,lab_size=6)
dev.off()
data_use = melt(res[,-c(1,2)])
pdf("miRNA_differ_R_rank_violinPlot.pdf")
ggplot(data_use, aes(x = variable, y = value))+ geom_violin(aes(fill = variable),scale = "width", alpha=0.7)+geom_boxplot(width = 0.1)+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), text = element_text(size = 12),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 20, face="bold", hjust=0.5))+ scale_y_continuous(name = "Pathway Rank")+scale_x_discrete(name = NULL)+guides(fill=FALSE)+ scale_fill_brewer(palette = "Set1")+stat_compare_means(label = "p.signif", method = "t.test", ref.group = "R0.2")
dev.off()


# lncRNA

load("lncRNA2Disease_all.RData")##	data_use
res = data_use
for(r in seq(0.1,0.9,0.1)){
print(r)
for(lnc in intersect(unique(res$ncRNA.Ensembl), lncRNA)){
load(paste0("/XXX/RWR", r, "/res", "/lncRNA_GSEA_res.RData"))
input = lncRNA_GSEA_res[[lnc]]
input[is.na(input$NES), "NOM p-val"] = 1
input[is.na(input$NES), "NES"] = -1
input[input$NES < 0, "NOM p-val"] = input[input$NES < 0, "NOM p-val"]+1
input.order = input[order(input$`NOM p-val`, -input$NES), ]
res[res$ncRNA.Ensembl == lnc, paste0("R", r)] = match(res[res$ncRNA.Ensembl == lnc, "Pathway"], input.order$NAME)
}
save(res, file = "lncRNA_differ_R_rank_comparison.RData")
}
res = na.omit(res)
cor(res[,-c(1,2)])
library(ggcorrplot)
corr <- round(cor(res[,-c(1,2)]), 1)
pdf("lncRNA_differ_R_rank_cor.pdf")
ggcorrplot(corr,type = "upper",outline.color = "white",lab=TRUE,lab_size=6)
dev.off()
data_use = melt(res[,-c(1,2)])
pdf("lncRNA_differ_R_rank_violinPlot.pdf")
ggplot(data_use, aes(x = variable, y = value))+ geom_violin(aes(fill = variable),scale = "width", alpha=0.7)+geom_boxplot(width = 0.1)+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), text = element_text(size = 12),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 20, face="bold", hjust=0.5))+ scale_y_continuous(name = "Target Pathway Rank")+scale_x_discrete(name = NULL)+guides(fill=FALSE)+ scale_fill_brewer(palette = "Set1")+stat_compare_means(label = "p.signif", method = "t.test", ref.group = "R0.2")
dev.off()


# circRNA

load("circRNA2DiseasePathway.RData")##	data_use
res = data_use
for(r in seq(0.1,0.9,0.1)){
print(r)
for(circ in intersect(unique(res$circRNA), circRNA)){
load(paste0("/XXX/RWR", r, "/res", "/circRNA_GSEA_res.RData"))
input = circRNA_GSEA_res[[circ]]
input[is.na(input$NES), "NOM p-val"] = 1
input[is.na(input$NES), "NES"] = -1
input[input$NES < 0, "NOM p-val"] = input[input$NES < 0, "NOM p-val"]+1
input.order = input[order(input$`NOM p-val`, -input$NES), ]
res[res$circRNA == circ, paste0("R", r)] = match(res[res$circRNA == circ, "Pathway"], input.order$NAME)
}
save(res, file = "circRNA_differ_R_rank_comparison.RData")
}
res = na.omit(res)
cor(res[,-c(1,2)])
library(ggcorrplot)
corr <- round(cor(res[,-c(1,2)]), 1)
pdf("circRNA_differ_R_rank_cor.pdf")
ggcorrplot(corr,type = "upper",outline.color = "white",lab=TRUE,lab_size=6)
dev.off()
data_use = melt(res[,-c(1,2)])
pdf("circRNA_differ_R_rank_violinPlot.pdf")
ggplot(data_use, aes(x = variable, y = value))+ geom_violin(aes(fill = variable),scale = "width", alpha=0.7)+geom_boxplot(width = 0.1)+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), text = element_text(size = 12),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 20, face="bold", hjust=0.5))+ scale_y_continuous(name = "Target Pathway Rank")+scale_x_discrete(name = NULL)+guides(fill=FALSE)+ scale_fill_brewer(palette = "Set1")+stat_compare_means(label = "p.signif", method = "t.test", ref.group = "R0.2")
dev.off()

load("miRNA_differ_R_rank_comparison.RData")
res = na.omit(res)
data_use = melt(res[,-c(1,2)])
data_use$ncRNA = "miRNA"
data_plot = data_use
load("lncRNA_differ_R_rank_comparison.RData")
res = na.omit(res)
data_use = melt(res[,-c(1,2)])
data_use$ncRNA = "lncRNA"
data_plot = rbind(data_plot, data_use)
load("circRNA_differ_R_rank_comparison.RData")
res = na.omit(res)
data_use = melt(res[,-c(1,2)])
data_use$ncRNA = "circRNA"
data_plot = rbind(data_plot, data_use)
data_plot$ncRNA = factor(data_plot$ncRNA, levels = unique(data_plot$ncRNA))
library(limma)
r = as.character(as.numeric(strsplit2(data_plot$variable, "[.]")[,2])*(0.1))
data_plot$variable = r
pdf("ncRNA_differ_R_rank_column_violinPlot.pdf", height = 12)
ggplot(data_plot, aes(x = variable, y = value))+ geom_violin(aes(fill = variable),scale = "width", alpha=0.7) +facet_wrap(~ncRNA, ncol=1, scales = "free_y") +geom_boxplot(width = 0.1)+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), text = element_text(size = 16, face = "bold"),axis.title = element_text(face="bold"), plot.title = element_text(size = 20, face="bold", hjust=0.5))+ scale_y_continuous(name = "Target Pathway Rank")+scale_x_discrete(name = "Restart Possibility")+guides(fill=FALSE)+scale_fill_brewer(palette = "Set1")+stat_compare_means(label = "p.signif", method = "t.test", ref.group = "0.2")
dev.off()
pdf("ncRNA_differ_R_rank_row_violinPlot.pdf", width = 12)
ggplot(data_plot, aes(x = variable, y = value))+ geom_violin(aes(fill = variable),scale = "width", alpha=0.7) +facet_wrap(~ncRNA, scales = "free_y") +geom_boxplot(width = 0.1)+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), text = element_text(size = 16, face = "bold"),axis.title = element_text(face="bold"), plot.title = element_text(size = 20, face="bold", hjust=0.5))+ scale_y_continuous(name = "Target Pathway Rank")+scale_x_discrete(name = "Restart Possibility")+guides(fill=FALSE)+scale_fill_brewer(palette = "Set1")+stat_compare_means(label = "p.signif", method = "t.test", ref.group = "0.2")
dev.off()
