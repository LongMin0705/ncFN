setwd("/XXX")
library(miRBaseConverter)
library(openxlsx)
library(dplyr)
## ncbi
NCBI_ID = read.table("Homo_sapiens.gene_info.txt",header = T,stringsAsFactors = F,fill = T, sep = "\t",quote = "")
NCBI_ID_sub = NCBI_ID[NCBI_ID$tax_id=="9606",c("GeneID", "Symbol", "dbXrefs")]
ensembl = c()
ensembl = apply(NCBI_ID_sub, 1, function(x){
if(grepl("Ensembl", x[3])){
temp = unlist(strsplit(as.character(x[3]), "[|]"))
ensembl = c(ensembl, unlist(strsplit(temp[grep("Ensembl", temp)], ":"))[2])
}else {ensembl = c(ensembl, NA)}
})
NCBI_ID_sub = unique(cbind(NCBI_ID_sub[,1:2], Ensembl=ensembl))

##  miRBase
miRNA_all <- read.xlsx("miRNA.xlsx")
temp <- miRNA_all[,c("Accession","ID","Mature1_Acc","Mature1_ID")]
colnames(temp) <- c("Accession","ID","Mature2_Acc","Mature2_ID")
temp <- rbind(temp, miRNA_all[,c("Accession","ID","Mature2_Acc","Mature2_ID")])
colnames(temp) <- c("pri.ACC","pri.name","ACC","name")
miRNA_name <- temp[grep("hsa",temp$name),]
miRNA_name <- na.omit(unique(miRNA_name))

##  HGNC
miRNA_HGNC = read.table("RNA_micro_HGNC.txt",header = T,stringsAsFactors = F, fill=T)
miRNA_HGNC_sub = miRNA_HGNC[,c("symbol","prev_name")]
multi = grep("\\|", miRNA_HGNC_sub$prev_name)
miRNA_HGNC_multi = data.frame()
for(i in multi){
  pre_name = unlist(strsplit(miRNA_HGNC_sub$prev_name[i], "\\|"))
  temp = data.frame(symbol=rep(miRNA_HGNC_sub$symbol[i],length(pre_name)),
                    prev_name=pre_name)
  miRNA_HGNC_multi = rbind(miRNA_HGNC_multi, temp)
}
miRNA_HGNC_sub2 = miRNA_HGNC_sub[-multi,]
miRNA_HGNC_sub2 = rbind(miRNA_HGNC_sub2, miRNA_HGNC_multi)
colnames(miRNA_HGNC_sub2)[2] = "pri.name"
miRNA_HGNC_use = na.omit(unique(merge(miRNA_name, miRNA_HGNC_sub2, by="pri.name", all.x=T)))

##  circBase
circRNA_ID <- read.table("hsa_hg19_circRNA.txt", header = T, sep = "\t", stringsAsFactors = F)
circRNA_info <- circRNA_ID[, c("circRNA.ID", "best.transcript", "gene.symbol")]
colnames(circRNA_info) <- c("circRNA_ID", "best_transcript", "gene_symbol")

layer_first <- read.table("big_network_nodes_symbol.txt", header = F, stringsAsFactors = F)
colnames(layer_first) <- c("entrez", "symbol", "entrez_1")

### 1 pathway interaction 375205 11555 #####
path_interact <- read.table("path_interact_new.txt", header = F, stringsAsFactors = F)
data <- unique(path_interact)
colnames(data) = c("EntrezA","EntrezB")
ID_list_1 <-NCBI_ID_sub[match(data[,1],NCBI_ID_sub$GeneID),1:2]
ID_list_2 <-NCBI_ID_sub[match(data[,2],NCBI_ID_sub$GeneID),1:2]
data = cbind(ID_list_1, ID_list_2)
colnames(data) = c("EntrezA","SymbolA","EntrezB","SymbolB")
path_interact <- data
delete <- apply(path_interact[, c(1,3)], 1, function(x){
  return(x[1] == x[2])
})##self-loop

### 2 protein-protein interaction 83065 12512 #####
ppi <- read.table("ppi_selected.txt", header = T, sep = "\t", stringsAsFactors = F)##  83065
temp = c()
for(i in seq(nrow(ppi))){
  temp[i] = 10-length(which(ppi[i,]==""))
}
sources = ppi[,-c(1,2)]
sources[which(sources == "", arr.ind=T)] = NA
sourcesl = c()
for(i in 1:ncol(sources)){
sourcesl = c(sourcesl, c(sources[,i]))
}
data <- as.matrix(unique(ppi[,1:2]))
delete <- apply(data, 1, function(x){
  return(x[1] == x[2])
})
ID_list_1 <- NCBI_ID_sub[match(data[,1],NCBI_ID_sub$GeneID),1:2]
ID_list_2 <- NCBI_ID_sub[match(data[,2],NCBI_ID_sub$GeneID),1:2]##  19 NA
ID_list_2[,1] <- data[,2]
RNA_NA = which(is.na(ID_list_2$Symbol))
ID_list_2[RNA_NA,2] = layer_first[match(ID_list_2[RNA_NA,1], layer_first$entrez),2]
data = cbind(ID_list_1, ID_list_2)
colnames(data) = c("EntrezA","SymbolA","EntrezB","SymbolB")
ppi <- data

### 3 Transcription factors-Target genes 4692 1961 #####
TF_gene <- read.table("transfac.txt", header = F, sep = "\t", stringsAsFactors = F)
data = TF_gene
delete <- apply(data, 1, function(x){
  return(x[1] == x[2])
})
data <- unique(data[!delete, ])
ID_list_1 <- NCBI_ID_sub[match(data[,1],NCBI_ID_sub$GeneID),1:2]
ID_list_2 <- NCBI_ID_sub[match(data[,2],NCBI_ID_sub$GeneID),1:2]
ID_list_2[,1] <- data[,2]
RNA_NA = which(is.na(ID_list_2$Symbol))
ID_list_2[RNA_NA,2] = layer_first[match(ID_list_2[RNA_NA,1], layer_first$entrez),2]
data = cbind(ID_list_1, ID_list_2)
colnames(data) = c("EntrezA","SymbolA","EntrezB","SymbolB")
TF_gene <- data

### 4 starbase lncRNA-mRNA, selection based on pathway interaction 38663 3440 10267 #####
lnc_RNA <- read.table("lncRNA-RNA.txt", header = T, sep = "\t", stringsAsFactors = F)
lnc_mRNA <- lnc_RNA[lnc_RNA$pairGeneType == "protein_coding", c(1,2,4,5)]
data <- unique(lnc_mRNA)
data$pairGeneEntrez = NA
for(i in 1:nrow(data)){
  id = NCBI_ID_sub[match(data[i, 1], NCBI_ID_sub$Ensembl), 2]
  if(!is.na(id)){
    if(id != data[i, 2]) data[i, 2] = id
  }
  id = NCBI_ID_sub[match(data[i, 2], NCBI_ID_sub$Symbol), 3]
  if(!is.na(id)){
    if(id != data[i, 1]) data[i, 1] = id
  }
  id = NCBI_ID_sub[match(data[i, 3], NCBI_ID_sub$Ensembl), 1]
  if(!is.na(id)) data[i, 5] = id
  id = NCBI_ID_sub[match(data[i, 4], NCBI_ID_sub$Symbol), 1]
  if(!is.na(id) & is.na(data[i, 5])) data[i, 5] = id
}
data1 <- unique(na.omit(data[,-3]))
RNA_ID = unique(data1[,1:2])
RNA_multi = names(which(table(RNA_ID$geneName)>1))
for(i in RNA_multi){
  data1 = data1[-which(data1[,2]==i),]
}
RNA_ID = unique(data1[,3:4])
RNA_multi = names(which(table(RNA_ID$pairGeneName)>1))
for(i in RNA_multi){
  data1 = data1[-which(data1[,3]==i),]
}
RNA_ID = unique(data1[,3:4])
RNA_multi = names(which(table(RNA_ID$pairGeneEntrez)>1))
for(i in RNA_multi){
  data1 = data1[-which(data1[,4]==i),]
}
lnc_mRNA_in <- unique(na.omit(data1))
lnc_mRNA_in[which(lnc_mRNA_in$geneName=="IGHA2"), 1] = "ENSG00000211890"
lnc_mRNA_in = lnc_mRNA_in[grep("ENSG", lnc_mRNA_in$geneID),]

### 5 LncRNA2Target lncRNA-mRNA 686 142 400 #####
l2t <- read.xlsx("lncRNA_target_from_low_throughput_experiments_V3.xlsx")
l2t_sub <- l2t[l2t$Species==9606,]
l2t_mi <- unique(l2t_sub[grep("MIR",l2t_sub$Target_official_symbol),c("Ensembl_ID","LncRNA_official_symbol","Target_official_symbol")])
l2t_m <- unique(l2t_sub[-grep("MIR",l2t_sub$Target_official_symbol),c("Ensembl_ID","LncRNA_official_symbol","Target_entrez_gene_ID","Target_official_symbol")])
data <- l2t_m
data = data[-which(is.na(data[,1]) & is.na(data[,2])), ]
data = data[grep("ENSG",data[,1]),]
data = data[-which(data[,3]==0),]
for(i in 1:nrow(data)){
  id = NCBI_ID_sub[match(data[i, 1], NCBI_ID_sub$Ensembl), 2]
  if(!is.na(id)){
    if(id != data[i, 2]) data[i, 2] = id
  }
  id = NCBI_ID_sub[match(data[i, 2], NCBI_ID_sub$Symbol), 3]
  if(!is.na(id)){
    if(id != data[i, 1]) data[i, 1] = id
  }
  id = NCBI_ID_sub[match(data[i, 3], NCBI_ID_sub$GeneID), 2]
  if(!is.na(id)){
    if(id != data[i, 4]) data[i, 4] = id
  }
  id = NCBI_ID_sub[match(data[i, 4], NCBI_ID_sub$Symbol), 1]
  if(!is.na(id)){
    if(id != data[i, 3]) data[i, 3] = id
  }
}
l2t_m_in = data

### 6 starbase  miRNA-mRNA 82073 613 8882 #####
miRNA_mRNA <- read.table("miRNA-mRNA.txt", header = T, stringsAsFactors = F)
data <- miRNA_mRNA
data <- data[data$degraExpNum>0&data$clipExpNum>0,1:4]
data$geneEntrez = NA
for(i in 1:nrow(data)){
  id = miRNA_name[match(data[i, 1], miRNA_name$ACC), 4]
  if(!is.na(id)){
    if(id != data[i, 2]) data[i, 2] = id
  }
  id = miRNA_name[match(data[i, 2], miRNA_name$name), 3]
  if(!is.na(id)){
    if(id != data[i, 1]) data[i, 1] = id
  }
  id = NCBI_ID_sub[match(data[i, 3], NCBI_ID_sub$Ensembl), 1]
  if(!is.na(id)) data[i, 5] = id
  id = NCBI_ID_sub[match(data[i, 4], NCBI_ID_sub$Symbol), 1]
  if(!is.na(id) & is.na(data[i, 5])) data[i, 5] = id
}
data1 = unique(na.omit(data[,-3]))
RNA_ID = unique(data1[,3:4])
RNA_multi = names(which(table(RNA_ID$geneName)>1))
for(i in RNA_multi){
  data1 = data1[-which(data1[,3]==i),]
}
RNA_ID = unique(data1[,3:4])
RNA_multi = names(which(table(RNA_ID$geneEntrez)>1))
for(i in RNA_multi){
  data1 = data1[-which(data1[,4]==i),]
}
miRNA_mRNA_1 <- unique(na.omit(data1))

### 7 mirTarBase miRNA-mRNA 8659 739 2848 ######
mirTarBase <- read.xlsx("hsa_MTI_mirTarBase.xlsx")
mirTarBase_sub <- unique(mirTarBase[-grep("Weak",mirTarBase$Support.Type),c("miRNA","Target.Gene.(Entrez.Gene.ID)","Target.Gene")])
data = mirTarBase_sub
ID_list_1 <- miRNA_name[match(data[,1], miRNA_name$name),3]
data <- unique(cbind(ID_list_1,data,stringsAsFactors=F))
RNA_NA = unique(data[which(is.na(data[,1])),])
RNA_NA[which(RNA_NA$miRNA=="hsa-miR-128b"),2] = "hsa-mir-128-2"
data1 = data.frame()
data2 = c()
for(i in 1:nrow(RNA_NA)){
  temp = miRNA_name[which(miRNA_name$pri.name==tolower(RNA_NA[i,2])),3:4]
  if(dim(temp)[1] != 0){
    temp = cbind(temp, RNA_NA[i,3:4])
    data1 = rbind(data1, temp)
  }else data2 = c(data2, RNA_NA[i,2])
}
RNA_NA[which(RNA_NA$miRNA=="hsa-miR-1254"),1] = "MIMAT0005905"
colnames(RNA_NA) = colnames(data1)
data1 = rbind(data1, RNA_NA[which(RNA_NA$name=="hsa-miR-1254"),])
data = na.omit(data)
colnames(data) = colnames(data1)
data = na.omit(unique(rbind(data,data1)))
for(i in 1:nrow(data)){
  id = miRNA_name[match(data[i, 1], miRNA_name$ACC), 4]
  if(!is.na(id)){
    if(id != data[i, 2]) data[i, 2] = id
  }
  id = miRNA_name[match(data[i, 2], miRNA_name$name), 3]
  if(!is.na(id)){
    if(id != data[i, 1]) data[i, 1] = id
  }
  id = NCBI_ID_sub[match(data[i, 3], NCBI_ID_sub$GeneID), 2]
  if(!is.na(id)){
    if(id != data[i, 4]) data[i, 4] = id
  }
  id = NCBI_ID_sub[match(data[i, 4], NCBI_ID_sub$Symbol), 1]
  if(!is.na(id)){
    if(id != data[i, 3]) data[i, 3] = id
  }
}
mirTarBase_in = na.omit(unique(data))
colnames(mirTarBase_in) = c("miRNAID","miRNAname","geneEntrez","geneSymbol")

### 8 TransMir TF-miRNA 3667  324(TF) 578#####
transMir <- read.csv("hsa_tsv_TransMir/hsa.tsv", sep = "\t",header = F, stringsAsFactors = F)
data <- na.omit(unique(transMir[transMir[,7] == "literature",1:2]))
ID_list_1 <- NCBI_ID_sub[match(data[,1],NCBI_ID_sub$Symbol),1:2]
ID_list_1 = cbind(ID_list_1, pri.name=data[,2])
ID_list_1 = na.omit(unique(ID_list_1))
ID_list_2 <- miRNA_name[miRNA_name$pri.name %in% unique(ID_list_1[,3]),]
RNA_NA = setdiff(unique(ID_list_1[,3]), ID_list_2$pri.name)
RNA_NA_dat = data.frame()
for(i in RNA_NA){
  versions = checkMiRNAVersion(i, verbose = F)
  temp = miRNA_NameToAccession(i,version = versions)
  colnames(temp) = c("miRNAname","ACC")
  RNA_NA_dat = rbind(RNA_NA_dat, temp)
}
RNA_NA_dat = na.omit(RNA_NA_dat)
RNA_NA_dat1 = data.frame()
RNA_NA_dat2 = data.frame()
for(i in seq(nrow(RNA_NA_dat))){
  if(grepl("MIMAT",RNA_NA_dat[i,2])){
    if(length(which(miRNA_name$ACC==RNA_NA_dat[i,2]))!=0){
      temp = miRNA_name[miRNA_name$ACC==RNA_NA_dat[i,2], ]
      temp$origial = RNA_NA_dat[i,1]
      RNA_NA_dat1 = rbind(RNA_NA_dat1,temp)
    }
  }
  else{
    if(length(which(miRNA_name$pri.ACC==RNA_NA_dat[i,2]))!=0){
      temp = miRNA_name[miRNA_name$pri.ACC==RNA_NA_dat[i,2], ]
      temp$origial = RNA_NA_dat[i,1]
      RNA_NA_dat2 = rbind(RNA_NA_dat2, temp)
    }
  }
}
RNA_NA_dat1 = unique(rbind(RNA_NA_dat1[,c(5,3,4)],RNA_NA_dat2[,c(5,3,4)]))
colnames(RNA_NA_dat1) = colnames(ID_list_2)[2:4]
ID_list_2 = unique(na.omit(rbind(ID_list_2[,2:4], RNA_NA_dat1)))
ID_list = merge(ID_list_1,ID_list_2,by.x="pri.name",all.x=T)
ID_list = na.omit(unique(ID_list[,2:5]))
data1 <- ID_list
transMir_in = unique(na.omit(data1))

### 9 LncBase lncRNA-miRNA 71 59(miRNA) 18(lncRNA) ####
lncBase <- read.table("LncBasev2_download.csv", sep = "\t",quote = "",header = T,fill = T,stringsAsFactors = F)
lncBase_sub <- lncBase[lncBase$species == "Homo sapiens" | lncBase$species == "Homo Sapiens", ]
lncBase_sub <- lncBase_sub[lncBase_sub$method == "Biotin-qPCR" |lncBase_sub$method == "Luciferase Reporter Assay" |lncBase_sub$method == "Northern Blot" | lncBase_sub$method == "AGO-IP" | lncBase_sub$method == "qPCR",]
data <- unique(lncBase_sub[,1:3])
data[which(data[,2]=="CCAT2"),1] = "ENSG00000280997"
data[which(data[,2]=="H19"),1] = "ENSG00000130600"
data[which(is.na(data[,2])),2] = "SNHG16"
data = unique(data[grep("ENSG",data[,1]),])
RNA_ID = unique(data[,1:2])
RNA_multi = names(which(table(RNA_ID$geneId)>1))
for(i in RNA_multi){
  data = data[-which(data[,1]==i),]
}
data$miRNAid <- miRNA_name[match(data$mirna, miRNA_name$name),3]
RNA_NA = rbind(data[which(is.na(data[,4])),1:2],data[which(is.na(data[,4])),1:2])
RNA_NA = cbind(RNA_NA, mirna = c("hsa-miR-217-5p","hsa-miR-217-3p"))
RNA_NA = cbind(RNA_NA, miRNAid = c("MIMAT0000274","MIMAT0037308"))
data = na.omit(unique(data))
data = rbind(data, RNA_NA)
for(i in 1:nrow(data)){
  id = NCBI_ID_sub[match(data[i, 1], NCBI_ID_sub$Ensembl), 2]
  if(!is.na(id)){
    if(id != data[i, 2]) data[i, 2] = id
  }
  id = NCBI_ID_sub[match(data[i, 2], NCBI_ID_sub$Symbol), 3]
  if(!is.na(id)){
    if(id != data[i, 1]) data[i, 1] = id
  }
  id = miRNA_name[match(data[i, 3], miRNA_name$name), 3]
  if(!is.na(id)){
    if(id != data[i, 4]) data[i, 4] = id
  }
  id = miRNA_name[match(data[i, 4], miRNA_name$ACC), 4]
  if(!is.na(id)){
    if(id != data[i, 3]) data[i, 3] = id
  }
}
lncBase_in = unique(na.omit(data))

### 10 starbase miRNA-lncRNA 1879 553 383 ####
miRNA_lncRNA <- read.table("miRNA-lncRNA.txt", header = T, sep = "\t", stringsAsFactors = F, fill = T)
data <- miRNA_lncRNA[,c(1,2,3,4,10,11)]
data <- data[data$degraExpNum>0 & data$clipExpNum>0,]
data = data[, 1:4]
for(i in 1:nrow(data)){
  id = miRNA_name[match(data[i, 1], miRNA_name$ACC), 4]
  if(!is.na(id)){
    if(id != data[i, 2]) data[i, 2] = id
  }
  id = miRNA_name[match(data[i, 2], miRNA_name$name), 3]
  if(!is.na(id)){
    if(id != data[i, 1]) data[i, 1] = id
  }
  id = NCBI_ID_sub[match(data[i, 3], NCBI_ID_sub$Ensembl), 2]
  if(!is.na(id)){
    if(id != data[i, 4]) data[i, 4] = id
  }
  id = NCBI_ID_sub[match(data[i, 4], NCBI_ID_sub$Symbol), 3]
  if(!is.na(id)){
    if(id != data[i, 3]) data[i, 3] = id
  }
}
miRNA_lncRNA_1 <- na.omit(unique(data))## 1879

### 11 starbase miRNA-circRNA 64880 642 7958 #####
miRNA_circRNA <- read.table("miRNA-circRNA.txt", header = T, sep = "\t", stringsAsFactors = F)
data <- miRNA_circRNA
data <- data[data$degraExpNum>0 & data$clipExpNum>0, 1:4]
ID_list_2 <- circRNA_info[match(data[,3], circRNA_info$best_transcript), ]
ID_list_2[is.na(ID_list_2$circRNA_ID), 1] <- data[grep("hsa", data[,3]), 3]
ID_list <- cbind(data[,1:2], circRNA_ID = ID_list_2[,1])
data <- unique(na.omit(ID_list))
for(i in 1:nrow(data)){
  id = miRNA_name[match(data[i, 1], miRNA_name$ACC), 4]
  if(!is.na(id)){
    if(id != data[i, 2]) data[i, 2] = id
  }
  id = miRNA_name[match(data[i, 2], miRNA_name$name), 3]
  if(!is.na(id)){
    if(id != data[i, 1]) data[i, 1] = id
  }
}
miRNA_circRNA_1 <- na.omit(unique(data))

### 12 lncRNA2target lncRNA-miRNA  1156 160(lnc) 443 #####
data = l2t_mi
data = data[grep("ENSG",data[,1]),]
for(i in 1:nrow(data)){
  id = NCBI_ID_sub[match(data[i, 1], NCBI_ID_sub$Ensembl), 2]
  if(!is.na(id)){
    if(id != data[i, 2]) data[i, 2] = id
  }
  id = NCBI_ID_sub[match(data[i, 2], NCBI_ID_sub$Symbol), 3]
  if(!is.na(id)){
    if(id != data[i, 1]) data[i, 1] = id
  }
}
ID_list_2 = data.frame()
for(i in data[,3]){
  temp = miRNA_HGNC_use[miRNA_HGNC_use$symbol==i,3:5]
  ID_list_2 = rbind(ID_list_2,temp)
}
colnames(ID_list_2) = c("miRNAid","miRNAname","Target_official_symbol")
ID_list = na.omit(unique(merge(data,ID_list_2,by="Target_official_symbol",all=T)))
l2t_mi_in <- unique(na.omit(ID_list[,-1]))

### ID 1 to 1 ####
##  lncRNA  Ensembl symbol
head(lnc_mRNA_in)
lncRNA_dat = lnc_mRNA_in[,1:2]
colnames(lncRNA_dat) = c("Ensembl","Symbol")
temp = l2t_m_in[,1:2]
colnames(temp) = c("Ensembl","Symbol")
lncRNA_dat = rbind(lncRNA_dat, temp)
temp = lncBase_in[,1:2]
colnames(temp) = c("Ensembl","Symbol")
lncRNA_dat = rbind(lncRNA_dat, temp)
temp = miRNA_lncRNA_1[,3:4]
colnames(temp) = c("Ensembl","Symbol")
lncRNA_dat = rbind(lncRNA_dat, temp)
temp = l2t_mi_in[,1:2]
colnames(temp) = c("Ensembl","Symbol")
lncRNA_dat = rbind(lncRNA_dat, temp)
lncRNA_dat = unique(lncRNA_dat)## 3585

names(which(table(lncRNA_dat$Ensembl)>1))##  2
NCBI_ID_sub[which(NCBI_ID_sub$Ensembl == "ENSG00000114374"),]
lncRNA_dat[which(lncRNA_dat$Ensembl=="ENSG00000114374"), 2] = "TTTY15"
head(l2t_mi_in)
lnc_mRNA_in[which(lnc_mRNA_in$geneID == "ENSG00000114374"),2] = "TTTY15"
l2t_m_in[which(l2t_m_in$Ensembl_ID == "ENSG00000114374"),2] = "TTTY15"
lncBase_in[which(lncBase_in$geneId == "ENSG00000114374"),2] = "TTTY15"
miRNA_lncRNA_1[which(miRNA_lncRNA_1$geneID == "ENSG00000114374"),4] = "TTTY15"
l2t_mi_in[which(l2t_mi_in$Ensembl_ID == "ENSG00000114374"),2] = "TTTY15"

NCBI_ID_sub[which(NCBI_ID_sub$Ensembl == "ENSG00000176124"),]
lncRNA_dat[which(lncRNA_dat$Ensembl=="ENSG00000176124"), 2] = "DLEU7-AS1"
lnc_mRNA_in[which(lnc_mRNA_in$geneID == "ENSG00000176124"),2] = "DLEU7-AS1"
l2t_m_in[which(l2t_m_in$Ensembl_ID == "ENSG00000176124"),2] = "DLEU7-AS1"
lncBase_in[which(lncBase_in$geneId == "ENSG00000176124"),2] = "DLEU7-AS1"
miRNA_lncRNA_1[which(miRNA_lncRNA_1$geneID == "ENSG00000176124"),4] = "DLEU7-AS1"
l2t_mi_in[which(l2t_mi_in$Ensembl_ID == "ENSG00000176124"),2] = "DLEU7-AS1"
lncRNA_dat = unique(lncRNA_dat)


##  miRNA ID Name
head(l2t_mi_in)
miRNA_dat = miRNA_mRNA_1[,1:2]
colnames(miRNA_dat) = c("ID","Name")
temp = mirTarBase_in[,1:2]
colnames(temp) = c("ID","Name")
miRNA_dat = rbind(miRNA_dat, temp)
temp = transMir_in[,3:4]
colnames(temp) = c("ID","Name")
miRNA_dat = rbind(miRNA_dat, temp)
temp = lncBase_in[,4:3]
colnames(temp) = c("ID","Name")
miRNA_dat = rbind(miRNA_dat, temp)
temp = miRNA_lncRNA_1[,1:2]
colnames(temp) = c("ID","Name")
miRNA_dat = rbind(miRNA_dat, temp)
temp = miRNA_circRNA_1[,1:2]
colnames(temp) = c("ID","Name")
miRNA_dat = rbind(miRNA_dat, temp)
temp = l2t_mi_in[,3:4]
colnames(temp) = c("ID","Name")
miRNA_dat = rbind(miRNA_dat, temp)
miRNA_dat = unique(miRNA_dat)

##  mRNA Entrez Symbol
head(l2t_m_in)
mRNA_dat = path_interact[,1:2]
colnames(mRNA_dat) = c("Entrez","Symbol")
temp = path_interact[,3:4]
colnames(temp) = c("Entrez","Symbol")
mRNA_dat = rbind(mRNA_dat, temp)
temp = ppi[,1:2]
colnames(temp) = c("Entrez","Symbol")
mRNA_dat = rbind(mRNA_dat, temp)
temp = ppi[,3:4]
colnames(temp) = c("Entrez","Symbol")
mRNA_dat = rbind(mRNA_dat, temp)
temp = TF_gene[,1:2]
colnames(temp) = c("Entrez","Symbol")
mRNA_dat = rbind(mRNA_dat, temp)
temp = TF_gene[,3:4]
colnames(temp) = c("Entrez","Symbol")
mRNA_dat = rbind(mRNA_dat, temp)
temp = lnc_mRNA_in[,4:3]
colnames(temp) = c("Entrez","Symbol")
mRNA_dat = rbind(mRNA_dat, temp)
temp = l2t_m_in[,3:4]
colnames(temp) = c("Entrez","Symbol")
mRNA_dat = rbind(mRNA_dat, temp)
temp = miRNA_mRNA_1[,4:3]
colnames(temp) = c("Entrez","Symbol")
mRNA_dat = rbind(mRNA_dat, temp)
temp = mirTarBase_in[,3:4]
colnames(temp) = c("Entrez","Symbol")
mRNA_dat = rbind(mRNA_dat, temp)
temp = transMir_in[,1:2]
colnames(temp) = c("Entrez","Symbol")
mRNA_dat = rbind(mRNA_dat, temp)
mRNA_dat = unique(mRNA_dat)

RNA_multi = names(which(table(mRNA_dat$Entrez)>1))
for(i in RNA_multi){
  mRNA_dat[which(mRNA_dat$Entrez == i),2] = NCBI_ID_sub[match(i, NCBI_ID_sub$GeneID),2]
}
RNA_ID = NCBI_ID_sub[match(RNA_multi, NCBI_ID_sub$GeneID),1:2]
for(i in 1:nrow(RNA_ID)){
  path_interact[which(path_interact$EntrezA==RNA_ID[i,1]),2] = RNA_ID[i,2]
  path_interact[which(path_interact$EntrezB==RNA_ID[i,1]),4] = RNA_ID[i,2]
  ppi[which(ppi$EntrezA==RNA_ID[i,1]),2] = RNA_ID[i,2]
  ppi[which(ppi$EntrezB==RNA_ID[i,1]),4] = RNA_ID[i,2]
  TF_gene[which(TF_gene$EntrezA==RNA_ID[i,1]),2] = RNA_ID[i,2]
  TF_gene[which(TF_gene$EntrezB==RNA_ID[i,1]),4] = RNA_ID[i,2]
  lnc_mRNA_in[which(lnc_mRNA_in$pairGeneEntrez==RNA_ID[i,1]),3] = RNA_ID[i,2]
  l2t_m_in[which(l2t_m_in$Target_entrez_gene_ID==RNA_ID[i,1]),4] = RNA_ID[i,2]
  miRNA_mRNA_1[which(miRNA_mRNA_1$geneEntrez==RNA_ID[i,1]),3] = RNA_ID[i,2]
  mirTarBase_in[which(mirTarBase_in$geneEntrez==RNA_ID[i,1]),4] = RNA_ID[i,2]
  transMir_in[which(transMir_in$GeneID==RNA_ID[i,1]),2] = RNA_ID[i,2]
}
mRNA_dat = unique(mRNA_dat)

RNA_multi = names(which(table(mRNA_dat$Symbol)>1))
RNA_ID = data.frame(Symbol = c("HSPA14","MEMO1"),GeneID = c("51182", "76890"), stringsAsFactors = F)
for(i in 1:nrow(RNA_ID)){
  mRNA_dat[which(mRNA_dat$Symbol==RNA_ID[i,1]),1] = RNA_ID[i,2]
  path_interact[which(path_interact$SymbolA==RNA_ID[i,1]),1] = RNA_ID[i,2]
  path_interact[which(path_interact$SymbolB==RNA_ID[i,1]),3] = RNA_ID[i,2]
  ppi[which(ppi$SymbolA==RNA_ID[i,1]),1] = RNA_ID[i,2]
  ppi[which(ppi$SymbolB==RNA_ID[i,1]),3] = RNA_ID[i,2]
  TF_gene[which(TF_gene$SymbolA==RNA_ID[i,1]),1] = RNA_ID[i,2]
  TF_gene[which(TF_gene$SymbolB==RNA_ID[i,1]),3] = RNA_ID[i,2]
  lnc_mRNA_in[which(lnc_mRNA_in$pairGeneName==RNA_ID[i,1]),4] = RNA_ID[i,2]
  l2t_m_in[which(l2t_m_in$Target_official_symbol==RNA_ID[i,1]),3] = RNA_ID[i,2]
  miRNA_mRNA_1[which(miRNA_mRNA_1$geneName==RNA_ID[i,1]),4] = RNA_ID[i,2]
  mirTarBase_in[which(mirTarBase_in$geneSymbol==RNA_ID[i,1]),3] = RNA_ID[i,2]
  transMir_in[which(transMir_in$Symbol==RNA_ID[i,1]),1] = RNA_ID[i,2]
}
mRNA_dat = unique(mRNA_dat)

##  circRNA 7958
circRNA_dat = unique(as.character(miRNA_circRNA_1$circRNA_ID))

write.table(path_interact, "Path_Interaction.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(ppi, "PPI.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(TF_gene, "TRANSFAC.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(lnc_mRNA_in, "StarBase_lncRNA_mRNA.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(l2t_m_in, "LncRNA2Target_lncRNA_mRNA.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(miRNA_mRNA_1, "StarBase_miRNA_mRNA.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(mirTarBase_in, "MirTarBase_miRNA_mRNA.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(transMir_in, "TransMir_mRNA_miRNA.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(lncBase_in, "LncBase_lncRNA_miRNA.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(miRNA_lncRNA_1, "StarBase_miRNA_lncRNA.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(miRNA_circRNA_1, "StarBase_miRNA_circRNA.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(l2t_mi_in, "LncRNA2Target_lncRNA_miRNA.txt", col.names = T, row.names = F, quote = F, sep = "\t")

### merge ####
data_net <- data.frame(ID_A = path_interact[,1], ID_B = path_interact[,3],stringsAsFactors = F)
temp <- data.frame(ID_A = ppi[,1], ID_B = ppi[,3],stringsAsFactors = F)
data_net <- unique(rbind(data_net,temp))
temp <- data.frame(ID_A = TF_gene[,1], ID_B = TF_gene[,3],stringsAsFactors = F)
data_net <- unique(rbind(data_net,temp))
temp <- data.frame(data1=paste(data_net[,1], data_net[,2]),data2=paste(data_net[,2], data_net[,1]))
net_rm <- match(temp[,1], temp[,2])
temp <- cbind(node=1:length(net_rm), net_rm)
uni_sel <- 1:dim(temp)[1]
for(i in 1:dim(temp)[1]){
  if((!is.na(temp[i,2])) & (temp[i,1]<temp[i,2])){
    uni_sel[temp[i,2]] <- temp[i,1]
  }
}
data_net_u <- data_net[unique(uni_sel),]
data_net <- na.omit(unique(data_net_u))
temp <- data.frame(ID_A = lnc_mRNA_in$geneID, ID_B = lnc_mRNA_in$pairGeneEntrez, stringsAsFactors = F)
data_net <- unique(rbind(data_net,temp))
length(na.omit(unique(union(data_net[,1],data_net[,2]))))
temp <- data.frame(ID_A = l2t_m_in$Ensembl_ID, ID_B = l2t_m_in$Target_entrez_gene_ID, stringsAsFactors = F)
data_net <- unique(rbind(data_net,temp))
length(na.omit(unique(union(data_net[,1],data_net[,2]))))
temp <- data.frame(ID_A = miRNA_mRNA_1$miRNAid, ID_B = miRNA_mRNA_1$geneEntrez, stringsAsFactors = F)
data_net <- unique(rbind(data_net,temp))
length(na.omit(unique(union(data_net[,1],data_net[,2]))))
temp <- data.frame(ID_A = mirTarBase_in$miRNAID, ID_B = mirTarBase_in$geneEntrez, stringsAsFactors = F)
data_net <- unique(rbind(data_net,temp))
length(na.omit(unique(union(data_net[,1],data_net[,2]))))
temp <- data.frame(ID_A = transMir_in$ACC, ID_B = transMir_in$GeneID, stringsAsFactors = F)
data_net <- unique(rbind(data_net,temp))
length(na.omit(unique(union(data_net[,1],data_net[,2]))))
temp <- data.frame(ID_A = lncBase_in$miRNAid, ID_B = lncBase_in$geneId, stringsAsFactors = F)
data_net <- unique(rbind(data_net,temp))
length(na.omit(unique(union(data_net[,1],data_net[,2]))))
temp <- data.frame(ID_A = miRNA_lncRNA_1$miRNAid, ID_B = miRNA_lncRNA_1$geneID, stringsAsFactors = F)
data_net <- unique(rbind(data_net,temp))
length(na.omit(unique(union(data_net[,1],data_net[,2]))))
temp <- data.frame(ID_A = miRNA_circRNA_1$miRNAid, ID_B = miRNA_circRNA_1$circRNA_ID, stringsAsFactors = F)
data_net <- unique(rbind(data_net,temp))
length(na.omit(unique(union(data_net[,1],data_net[,2]))))
temp <- data.frame(ID_A = l2t_mi_in$miRNAid, ID_B = l2t_mi_in$Ensembl_ID, stringsAsFactors = F)
data_net <- unique(rbind(data_net,temp))
length(na.omit(unique(union(data_net[,1],data_net[,2]))))
temp <- data.frame(data1=paste(data_net[,1], data_net[,2]),
                   data2=paste(data_net[,2], data_net[,1]))
net_rm <- match(temp[,1], temp[,2])
temp <- cbind(node=1:length(net_rm), net_rm)
uni_sel <- 1:dim(temp)[1]
for(i in 1:dim(temp)[1]){
  if((!is.na(temp[i,2])) & (temp[i,1]<temp[i,2])){
    uni_sel[temp[i,2]] <- temp[i,1]
  }
}
data_net_u <- data_net[unique(uni_sel),]
data_net_raw = data_net_u
data_net_raw_nodes = union(data_net_u[,1],data_net_u[,2])

### max component ####
library(igraph)
net.data = data_net_raw
colnames(net.data) = c("from", "to")
net <- graph_from_data_frame(d = net.data, directed = F)
all_component <- components(net)
nodes = names(V(net)[all_component[[1]]==1])
nodes_rm = setdiff(data_net_raw_nodes, nodes)
subnet = delete_vertices(net, nodes_rm)
all_component <- components(subnet)
data_net_max = as_edgelist(subnet)
data_net_max = as.data.frame(data_net_max, stringsAsFactors = F)
colnames(data_net_max) = c("ID_A", "ID_B")
data_net_max_nodes = union(data_net_max$ID_A,data_net_max$ID_B)
write.table(data_net_max, "network.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(data_net_max_nodes, "nodes.txt", sep = "\t", quote = F, col.names = F, row.names = F)

### interaction in network  #####
data_net_max_paste = paste(data_net_max[,1], data_net_max[,2], sep = " ")

path_interact.in = rep(0, length(data_net_max_paste))
res1 = paste(path_interact[,1], path_interact[,3], sep = " ")
res2 = paste(path_interact[,3], path_interact[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    path_interact.in[i] = 1
}
write.table(path_interact.in, "path_interact_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

ppi.in = rep(0, length(data_net_max_paste))
res1 = paste(ppi[,1], ppi[,3], sep = " ")
res2 = paste(ppi[,3], ppi[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    ppi.in[i] = 1
}
write.table(ppi.in, "ppi_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

dim(TF_gene)
TF_gene.in = rep(0, length(data_net_max_paste))
res1 = paste(TF_gene[,1], TF_gene[,3], sep = " ")
res2 = paste(TF_gene[,3], TF_gene[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    TF_gene.in[i] = 1
}
write.table(TF_gene.in, "TF_gene_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

SB_lnc_m.in = rep(0, length(data_net_max_paste))
res1 = paste(lnc_mRNA_in[,1], lnc_mRNA_in[,4], sep = " ")
res2 = paste(lnc_mRNA_in[,4], lnc_mRNA_in[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    SB_lnc_m.in[i] = 1
}
write.table(SB_lnc_m.in, "SB_lnc_m_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

L2T_lnc_m.in = rep(0, length(data_net_max_paste))
res1 = paste(l2t_m_in[,1], l2t_m_in[,3], sep = " ")
res2 = paste(l2t_m_in[,3], l2t_m_in[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    L2T_lnc_m.in[i] = 1
}
write.table(L2T_lnc_m.in, "L2T_lnc_m_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

SB_mi_m.in = rep(0, length(data_net_max_paste))
res1 = paste(miRNA_mRNA_1[,1], miRNA_mRNA_1[,4], sep = " ")
res2 = paste(miRNA_mRNA_1[,4], miRNA_mRNA_1[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    SB_mi_m.in[i] = 1
}
write.table(SB_mi_m.in, "SB_mi_m_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

MTB_mi_m.in = rep(0, length(data_net_max_paste))
res1 = paste(mirTarBase_in[,1], mirTarBase_in[,3], sep = " ")
res2 = paste(mirTarBase_in[,3], mirTarBase_in[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    MTB_mi_m.in[i] = 1
}
write.table(MTB_mi_m.in, "MTB_mi_m_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

TM_m_mi.in = rep(0, length(data_net_max_paste))
res1 = paste(transMir_in[,1], transMir_in[,3], sep = " ")
res2 = paste(transMir_in[,3], transMir_in[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    TM_m_mi.in[i] = 1
}
write.table(TM_m_mi.in, "TM_m_mi_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

LB_lnc_mi.in = rep(0, length(data_net_max_paste))
res1 = paste(lncBase_in[,1], lncBase_in[,4], sep = " ")
res2 = paste(lncBase_in[,4], lncBase_in[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    LB_lnc_mi.in[i] = 1
}
write.table(LB_lnc_mi.in, "LB_lnc_mi_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

SB_mi_lnc.in = rep(0, length(data_net_max_paste))
res1 = paste(miRNA_lncRNA_1[,1], miRNA_lncRNA_1[,3], sep = " ")
res2 = paste(miRNA_lncRNA_1[,3], miRNA_lncRNA_1[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    SB_mi_lnc.in[i] = 1
}
write.table(SB_mi_lnc.in, "SB_mi_lnc_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

SB_mi_circ.in = rep(0, length(data_net_max_paste))
res1 = paste(miRNA_circRNA_1[,1], miRNA_circRNA_1[,3], sep = " ")
res2 = paste(miRNA_circRNA_1[,3], miRNA_circRNA_1[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    SB_mi_circ.in[i] = 1
}
write.table(SB_mi_circ.in, "SB_mi_circ_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

temp = data_net_max[which(path_interact.in == 1),]
nodes = union(temp$ID_A,temp$ID_B)

L2T_lnc_mi.in = rep(0, length(data_net_max_paste))
res1 = paste(l2t_mi_in[,1], l2t_mi_in[,3], sep = " ")
res2 = paste(l2t_mi_in[,3], l2t_mi_in[,1], sep = " ")
for(i in 1:length(data_net_max_paste)){
  if(data_net_max_paste[i] %in% res1 | data_net_max_paste[i] %in% res2)
    L2T_lnc_mi.in[i] = 1
}
write.table(L2T_lnc_mi.in, "L2T_lnc_mi_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)

### info summary  ####
##  lncRNA
lncRNA <- data_net_max_nodes[grep("ENSG", data_net_max_nodes)]
lncRNA_index = match(lncRNA, data_net_max_nodes)-1
write.table(lncRNA, "lncRNA.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(lncRNA_index, "lncRNA_index.txt", sep = "\t", quote = F, col.names = F, row.names = F)
lncRNA_ID = lncRNA_dat[match(lncRNA, lncRNA_dat$Ensembl), ]
write.table(lncRNA_ID, "lncRNA_ID.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##  circRNA
circRNA <- data_net_max_nodes[grep("hsa_circ",data_net_max_nodes)]
circRNA_index = match(circRNA, data_net_max_nodes)-1
write.table(circRNA, "circRNA.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(circRNA_index, "circRNA_index.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(circRNA_info, "circRNA_ID.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##  miRNA
miRNA <- data_net_max_nodes[grep("MIMAT",data_net_max_nodes)]
miRNA_index = match(miRNA, data_net_max_nodes)-1
write.table(miRNA, "miRNA.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(miRNA_index, "miRNA_index.txt", sep = "\t", quote = F, col.names = F, row.names = F)
miRNA_ID = miRNA_dat[match(miRNA, miRNA_dat$ID), ]
write.table(miRNA_ID, "miRNA_ID.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##  mRNA
mRNA <- setdiff(data_net_max_nodes, c(lncRNA, circRNA, miRNA))
mRNA_index = match(mRNA, data_net_max_nodes)-1
write.table(mRNA, "mRNA.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(mRNA_index, "mRNA_index.txt", sep = "\t", quote = F, col.names = F, row.names = F)
mRNA_ID = mRNA_dat[match(mRNA, mRNA_dat$Entrez), ]
write.table(mRNA_ID, "mRNA_ID.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##  degree  ####
degree.all = as.matrix(table(c(data_net_max$ID_A, data_net_max$ID_B)))
write.table(degree.all, "degree.txt", sep = "\t", quote = F, col.names = FALSE, row.names = TRUE)
degree.mRNA = as.matrix(degree.all[mRNA,])
write.table(degree.mRNA, "degree_mRNA.txt", sep = "\t", quote = F, col.names = FALSE, row.names = TRUE)

data = table(c(data_net_max$ID_A, data_net_max$ID_B))
degree_stat = matrix(degree,nrow=length(degree))
all_degree_stat = cbind(as.numeric(names(degree)),as.numeric(degree_stat))
colnames(all_degree_stat) = c("degree","nodes.number")
all_degree_stat = as.data.frame(all_degree_stat)
reg = lm(log(all_degree_stat[,2])~log(all_degree_stat[,1]))
cozf = coef(reg)
power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
pdf("network_nodes_degree_distribution.pdf")
plot(all_degree_stat[,2]~all_degree_stat[,1],log="xy",xlab = "Node Degree", ylab = "Number of nodes", col = "#32CCFE",pch=16)
curve(power.law.fit,col="red",add=T,lty=2)
a = exp(cozf[[1]])
alpha = -cozf[[2]]
R.square = summary(reg)$r.squared
text_infor=paste("R square:",round(R.square,3),sep=" ")
text(500,1500,text_infor,cex =1.5)
dev.off()
print(paste("Alpha =", round(alpha, 3)))
print(paste("R square =", round(R.square, 3)))

data1 = data[mRNA]
degree = table(data1)
degree_stat = matrix(degree,nrow=length(degree))
all_degree_stat = cbind(as.numeric(names(degree)),as.numeric(degree_stat))
colnames(all_degree_stat) = c("degree","nodes.number")
all_degree_stat = as.data.frame(all_degree_stat)
reg = lm(log(all_degree_stat[,2])~log(all_degree_stat[,1]))
cozf = coef(reg)
power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
pdf("network_mRNA_degree_distribution.pdf")
plot(all_degree_stat[,2]~all_degree_stat[,1],log="xy",xlab = "mRNA Degree", 
     ylab = "Number of nodes", col = "#32CCFE",pch=16)#pch??¦Ì??
curve(power.law.fit,col="red",add=T,lty=2)
a = exp(cozf[[1]])
alpha = -cozf[[2]]
R.square = summary(reg)$r.squared
text_infor=paste("R square:",round(R.square,3),sep=" ")
text(400,600,text_infor,cex =1.5)
dev.off()

data1 = data[miRNA]
degree = table(data1)
degree_stat = matrix(degree,nrow=length(degree))
all_degree_stat = cbind(as.numeric(names(degree)),as.numeric(degree_stat))
colnames(all_degree_stat) = c("degree","nodes.number")
all_degree_stat = as.data.frame(all_degree_stat)
reg = lm(log(all_degree_stat[,2])~log(all_degree_stat[,1]))
cozf = coef(reg)
power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
pdf("network_miRNA_degree_distribution.pdf")
plot(all_degree_stat[,2]~all_degree_stat[,1],log="xy",xlab = "miRNA Degree", 
     ylab = "Number of nodes", col = "#32CCFE",pch=16)
curve(power.law.fit,col="red",add=T,lty=2)
a = exp(cozf[[1]])
alpha = -cozf[[2]]
R.square = summary(reg)$r.squared
text_infor=paste("R square:",round(R.square,3),sep=" ")
text(100,100,text_infor,cex =1.5)
dev.off()

data1 = data[lncRNA]
degree = table(data1)
degree_stat = matrix(degree,nrow=length(degree))
all_degree_stat = cbind(as.numeric(names(degree)),as.numeric(degree_stat))
colnames(all_degree_stat) = c("degree","nodes.number")
all_degree_stat = as.data.frame(all_degree_stat)
reg = lm(log(all_degree_stat[,2])~log(all_degree_stat[,1]))
cozf = coef(reg)
power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
pdf("network_lncRNA_degree_distribution.pdf")
plot(all_degree_stat[,2]~all_degree_stat[,1],log="xy",xlab = "lncRNA Degree", 
     ylab = "Number of nodes", col = "#32CCFE",pch=16)
curve(power.law.fit,col="red",add=T,lty=2)
a = exp(cozf[[1]])
alpha = -cozf[[2]]
R.square = summary(reg)$r.squared
text_infor=paste("R square:",round(R.square,3),sep=" ")
text(100,600,text_infor,cex =1.5)
dev.off()

data1 = data[circRNA]
degree = table(data1)
degree_stat = matrix(degree,nrow=length(degree))
all_degree_stat = cbind(as.numeric(names(degree)),as.numeric(degree_stat))
colnames(all_degree_stat) = c("degree","nodes.number")
all_degree_stat = as.data.frame(all_degree_stat)
reg = lm(log(all_degree_stat[,2])~log(all_degree_stat[,1]))
cozf = coef(reg)
power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
pdf("network_ciriRNA_degree_distribution.pdf")
plot(all_degree_stat[,2]~all_degree_stat[,1],log="xy",xlab = "ciriRNA Degree", 
     ylab = "Number of nodes", col = "#32CCFE",pch=16)
curve(power.law.fit,col="red",add=T,lty=2)
a = exp(cozf[[1]])
alpha = -cozf[[2]]
R.square = summary(reg)$r.squared
text_infor=paste("R square:",round(R.square,3),sep=" ")
text(50,600,text_infor,cex =1.5)
dev.off()
