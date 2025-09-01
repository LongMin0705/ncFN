# GSEA

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("/XXX")
lncRNA = read.table("lncRNA.txt", sep = "\t", header = FALSE)
mRNA = read.table("mRNA.txt", sep = "\t", header = FALSE)
miRNA = read.table("miRNA.txt", sep = "\t", header = FALSE)
circRNA = read.table("circRNA.txt", sep = "\t", header = FALSE)
ncRNA = union(miRNA[,1],union(lncRNA[,1],circRNA[,1]))

r = 0.2
run = paste0("/mnt/GSEA_Linux_4.0.2/gsea-cli.sh GSEAPreranked -gmx kegg_pathway.gmt -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk /XXX/RWR",r,"/rnk/replace_finalScore.rnk -scoring_scheme weighted -rpt_label replace -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 0 -rnd_seed timestamp -set_max 500 -set_min 1 -zip_report false -out /XXX/RWR",r,"/GSEA")
run_sh = data.frame(row.names = c(miRNA[,1],lncRNA[,1],circRNA[,1]))
for(i in c(miRNA[,1],lncRNA[,1],circRNA[,1])){
run_sh[i, "run"] = gsub("replace", i, run)
}
write.table(run_sh, "run.sh", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

run_sh = data.frame(row.names = c(miRNA[,1]))
for(i in c(miRNA[,1])){
run_sh[i, "run"] = gsub("replace", i, run)
}
write.table(run_sh, "run_mir.sh", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

run_sh = data.frame(row.names = c(lncRNA[,1]))
for(i in c(lncRNA[,1])){
run_sh[i, "run"] = gsub("replace", i, run)
}
write.table(run_sh, "run_lnc.sh", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

run_sh = data.frame(row.names = c(circRNA[,1]))
for(i in c(circRNA[,1])){
run_sh[i, "run"] = gsub("replace", i, run)
}
write.table(run_sh, "run_circ.sh", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# GSEA result summary
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

# change r
r = 0.2

# miRNA
path = paste0("/XXX/",r,"/GSEA")
setwd(path)
files = dir()[grep("^MIMAT", dir())]
miRNA_GSEA_res = list()
for(i in files){
setwd(paste(path, i, sep="/"))
name = unlist(strsplit(i, "[.]"))
pos = as.data.frame(fread(paste0("gsea_report_for_na_pos_", name[3], ".xls")))
pos.order = pos[order(pos$`NOM p-val`, -pos$NES), ]
neg = as.data.frame(fread(paste0("gsea_report_for_na_neg_", name[3], ".xls")))
if(neg[1,1] == "NAME"){
miRNA_GSEA_res[[name[1]]] = pos.order
}else{
neg.order = neg[order(neg$`NOM p-val`, -neg$NES), ]
miRNA_GSEA_res[[name[1]]] = rbind(pos.order, neg.order)
}
}
table(unlist(lapply(miRNA_GSEA_res, nrow)))
save(miRNA_GSEA_res, file = "miRNA_GSEA_res.RData")

# lncRNA
files = dir()[grep("^ENSG", dir())]
lncRNA_GSEA_res = list()
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
lncRNA_GSEA_res[[name[1]]] = input.order
}
table(unlist(lapply(lncRNA_GSEA_res, nrow)))
save(lncRNA_GSEA_res, file = "lncRNA_GSEA_res.RData")

# circRNA
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