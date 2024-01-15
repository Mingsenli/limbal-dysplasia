
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(dplyr)
library(tidyverse)
library(patchwork)
library(harmony)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(stringr)
library(data.table)
plan("multicore", workers = 48) ###set the compute core
options(future.globals.maxSize = 10000 * 1024^2)

library(cell_typeProfiler)
library(msigdbr)



dir = c('dysplasia/', 'adult1/', 'adult2/','adult3/','adult4/')

names(dir) = c('dysplasia','adult1', 'adult2','adult3', 'adult4') 


counts <- Read10X(data.dir =dir)
dysplasia_adult = CreateSeuratObject(counts,min.cells = 10, min.features = 300)

rm(counts,dir)

saveRDS(dysplasia_adult,file="./dysplasia_adult.SeuratObject.rds")

rm(dysplasia_adult)


dysplasia_adult<-readRDS(file="./dysplasia_adult.SeuratObject.rds")

table(dysplasia_adult@meta.data$orig.ident)

table(cornea_develop.anchor.sct.final@meta.data$orig.ident)

dysplasia_adult[["percent.mt"]] <- PercentageFeatureSet(dysplasia_adult, pattern = "^MT-")
hist(dysplasia_adult[["percent.mt"]]$percent.mt) 


HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(dysplasia_adult@assays$RNA)) 
HB.genes <- rownames(dysplasia_adult@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
dysplasia_adult[["percent.HB"]]<-PercentageFeatureSet(dysplasia_adult, features=HB.genes) 


VlnPlot(dysplasia_adult, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 2) #nFeature基因数中位数在1000左右，nCount在1万以下，主要在5/6千左右，


dysplasia_adult_filter <- subset(dysplasia_adult, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & percent.mt < 15 & percent.HB < 5)
dysplasia_adult@meta.data[1:5,]


dysplasia_adult_filter@meta.data[which(dysplasia_adult_filter@meta.data$orig.ident =="adult1"),"group"] <- "adult"
dysplasia_adult_filter@meta.data[which(dysplasia_adult_filter@meta.data$orig.ident =="adult2"),"group"] <- "adult"
dysplasia_adult_filter@meta.data[which(dysplasia_adult_filter@meta.data$orig.ident =="adult3"),"group"] <- "adult"
dysplasia_adult_filter@meta.data[which(dysplasia_adult_filter@meta.data$orig.ident =="adult4"),"group"] <- "adult"
dysplasia_adult_filter@meta.data[which(dysplasia_adult_filter@meta.data$orig.ident =="dysplasia"),"group"] <- "dysplasia"



saveRDS(dysplasia_adult_filter,file="./dysplasia_adult.filter.rds")


dysplasia_adult_filter <-readRDS(file="./dysplasia_adult.filter.rds")

dysplasia_adult.list <- SplitObject(dysplasia_adult_filter, split.by = "orig.ident")


for (i in 1:length(dysplasia_adult.list)) {
  
  dysplasia_adult.list[[i]] <- SCTransform(dysplasia_adult.list[[i]], vst.flavor = "v2", verbose = FALSE,variable.features.n = 3000,
                                           vars.to.regress = c("percent.mt","S.Score","G2M.Score")) }

features <- SelectIntegrationFeatures(object.list = dysplasia_adult.list, nfeatures = 3000)
dysplasia_adult.list <- PrepSCTIntegration(object.list = dysplasia_adult.list, anchor.features = features)

dysplasia_adult.anchor.sct <- FindIntegrationAnchors(object.list = dysplasia_adult.list, normalization.method = "SCT",
                                                     anchor.features = features,k.anchor = 10)

dysplasia_adult.anchor.sct <- IntegrateData(anchorset = dysplasia_adult.anchor.sct, normalization.method = "SCT")%>%
  RunPCA(verbose = FALSE,npcs=50)


saveRDS(dysplasia_adult.anchor.sct,file="./dysplasia_adult.anchor10.sct.rds")

dysplasia_adult.anchor.sct <-readRDS(file="./dysplasia_adult.anchor10.sct.rds")

ElbowPlot(dysplasia_adult.anchor.sct,ndims = 40,reduction="pca")

num=1:30

dysplasia_adult.anchor.sct.final <- dysplasia_adult.anchor.sct %>%
  
  RunUMAP(reduction = "pca", dims = num, verbose = FALSE) %>% 
  RunTSNE(reduction = "pca", dims = num, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = num) %>%
  FindClusters(resolution = 1.5,verbose = FALSE)


DefaultAssay(dysplasia_adult.anchor.sct.final) <- "RNA"

dysplasia_adult.anchor.sct.final <- NormalizeData(dysplasia_adult.anchor.sct.final, normalization.method = "LogNormalize", scale.factor = 10000)

DefaultAssay(dysplasia_adult.anchor.sct.final) <- "SCT"


saveRDS(dysplasia_adult.anchor.sct.final,file="./dysplasia_adult.anchor.sct.final.k.anchor10.rds")


DefaultAssay(dysplasia_adult.anchor.sct.final) <- "RNA"


cors = c("#20D3D3","#A0CBE8",'#99C945',"#B15928","#F28E2B","9900CC","#B6992D","#E15759","#79706E","#DDD822",
         "#D37295","#06918D",'#DAA51B',"#911958","#F0027F","#9D7660","#7A20EF","#E31A1C", "#153FB7",
         "#22A514","#345917", '#ED645A',"#223D6C",'#000000')


DimPlot(dysplasia_adult.anchor.sct.final,reduction = "umap",label = F,cols = cors)


dysplasia_adult.anchor.sct.final <-RenameIdents(dysplasia_adult.anchor.sct.final, "23"= "LSC","25"="LSC","7"="cornea.basal", 
                                                
                                                "6"="limbal.suprabasal","26"="limbal.suprabasal","9"="limbal.superfcial", 
                                                
                                                "15"="cornea.suprabasal","2"="cornea.suprabasal","29"="cornea.suprabasal",
                                                
                                                "12"="cornea.superfcial","10"="cornea.superfcial","21"="CESC",
                                                
                                                "0"="conjunctival.suprabasal", "4"="conjunctival.suprabasal","33"="conjunctival.suprabasal",
                                                
                                                "40"="conjunctival.suprabasal","13"="conjunctival.superfcial", "18"="conjunctival.superfcial",
                                                
                                                "19"="limbal.keratocytes", "14"="limbal.keratocytes","24"="limbal.keratocytes","20"="limbal.keratocytes",
                                                
                                                "30"="limbal.keratocytes","16"="limbal.fibroblasts","38"="limbal.fibroblasts","3"="cornea.keratocytes",
                                                
                                                "28"="cornea.keratocytes","22"="cornea.keratocytes","8"="cornea.keratocytes","41"="cornea.keratocytes",
                                                
                                                "32"="cornea.endothelium","11"="dysplasia1","44"="dysplasia1","1"="dysplasia2","17"="dysplasia3","5"="dysplasia4",
                                                
                                                "35"="melanocytes","31"= "macrophage","36"="immune.cell", "27"= "blood.vessels","37"= "blood.vessels", 
                                                
                                                "39"= "lymphatic.vessels","43"="unknown","45"="unknown","34"="unknown","42"="unknown")


dysplasia_adult.anchor.sct.final$cell_type<-Idents(dysplasia_adult.anchor.sct.final) 

Idents(dysplasia_adult.anchor.sct.final)=dysplasia_adult.anchor.sct.final@meta.data$cell_type

table(dysplasia_adult.anchor.sct.final$cell_type,dysplasia_adult.anchor.sct.final@meta.data$group)


all_genes <- rownames(dysplasia_adult.anchor.sct.final@assays$RNA@data)
keep_genes <- all_genes[!(grepl("^MT-", all_genes) | grepl("^RPL", all_genes) | grepl("^RPS", all_genes))]
dysplasia_adult.anchor.sct.final <- subset(dysplasia_adult.anchor.sct.final,features = keep_genes)



saveRDS(dysplasia_adult.anchor.sct.final,file="./dysplasia_adult.anchor.sct.final.k.anchor10.anno.normalized_new.rds")




###cell percentage

pB2_df <- table(dysplasia_adult.anchor.sct.final@meta.data$cell_type,dysplasia_adult.anchor.sct.final@meta.data$group) %>% melt()

colnames(pB2_df) <- c("cell_type","group","Number")

cell_type=c("LSC","limbal.suprabasal", "limbal.superfcial", 
            "cornea.basal","cornea.suprabasal","cornea.superfcial","cornea.endothelium",
            "CESC","conjunctival.suprabasal","conjunctival.superfcial",
            "limbal fbroblasts","limbal.keratocytes","cornea.keratocytes","melanocytes",
            "blood.vessels","lymphatic.vessels","macrophage","immune.cell",
            "dysplasia1","dysplasia2","dysplasia3","dysplasia4")


pB2_df$cell_type <- factor(pB2_df$cell_type,levels = cell_type)

pB2_df$group <- factor(pB2_df$group,levels = c("adult","dysplasia"))



pB2_df=na.omit( pB2_df)

pB4 <- ggplot(data = pB2_df, aes(x =Number, y =group, fill =  cell_type)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(y="",x="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB4

ggsave("dyplasia-percentage.pdf",pB4, width = 8, height = 7)

###############################




##differential gene expression analysis

type=c("LSC","cornea.basal","limbal.fibroblasts","limbal.keratocytes","cornea.keratocytes")


r.deg=data.frame()


for (i in 1:5) {
  
  deg=FindMarkers(dysplasia_adult.anchor.sct.final,ident.1 = "dysplasia",ident.2 = "adult",
                  logfc.threshold = 0.5, min.pct = 0.5,test.use = "wilcox",slot = "data",
                  group.by = "group",subset.ident =type[i],assay = "RNA")
  
  write.csv(deg,file = paste0( type[i],'_deg.csv') )
  deg$cell_type=type[i]
  deg$unm=i-1
  deg$gene=rownames(deg)
  r.deg=rbind(deg,r.deg)
  
}


for (i in 1:5) {
  
  diffgene <- subset(r.deg,subset=cell_type==type[i] & threshold=="Down")$gene
  
  
  erich.go.BP = enrichGO(gene =diffgene,
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = "BP",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
  
  write.csv(erich.go.BP, file = paste0(type[i],'_down.BP.csv'),sep = "\t",col.names = NA)
  
}


for (i in 1:5) {
  
  diffgene <- subset(r.deg,subset=cell_type==type[i] & threshold=="Up")$gene
  
  erich.go.BP = enrichGO(gene =diffgene,
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = "BP",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
  
  write.csv(erich.go.BP, file = paste0(type[i],'_up.BP.csv'),sep = "\t",col.names = NA)
  
}
#############################

##SCENIC analysis
pyscenic grn \
--num_workers 48 -m grnboost2 \
-o dysplasia_epithelium.loom.grn.output.tsv \
--seed 123 \
dysplasia_epithelium_filtered_data.loom \
/data/User/limingsen/data/RNA/singlecell/cisTarget_databases/human_all_TF.txt



pyscenic ctx \
dysplasia_epithelium.loom.grn.output.tsv \
/data/User/limingsen/data/RNA/singlecell/cisTarget_databases/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
/data/User/limingsen/data/RNA/singlecell/cisTarget_databases/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname /data/User/limingsen/data/RNA/singlecell/cisTarget_databases/motifs-v9nr_clust-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname dysplasia_epithelium_filtered_data.loom \
--mode "dask_multiprocessing" \
--output dysplasia_epithelium.loom.regulons.csv \
--num_workers 24 


srun -c 12 pyscenic aucell \
dysplasia_epithelium_all_data.loom \
dysplasia_epithelium.loom.regulons.csv \
-o dysplasia_epithelium.AUC.loom \
--num_workers 12 \
--seed 123

################################


####dynamo analysis in python
adata = sc.read_h5ad('dysplasia_adult_count_merge_dynamo_input.h5ad')

%matplotlib inline
from dynamo.dynamo_logger import main_info, LoggerManager
LoggerManager.main_logger.setLevel(LoggerManager.INFO)

dyn.pp.recipe_monocle(adata)
dyn.tl.dynamics(adata, cores=48)
dyn.tl.reduceDimension(adata)

dyn.tl.cell_velocities(adata,basis='pca')
dyn.vf.VectorField(adata,basis='pca', M=1000,cores=48)

dyn.tl.cell_velocities(adata,basis='umap')
dyn.vf.VectorField(adata,basis='umap', M=1000,cores=48)

