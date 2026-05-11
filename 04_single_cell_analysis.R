options(stringsAsFactors = F)
#加载R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(future)
library(future.apply)
plan(multisession, workers=4) ###set the compute core
options(future.globals.maxSize = 12000 * 1024^2)

#读取数据
#1、批量读取单细胞的数据
dir_name=c('Normal1','Normal2','Sepsis1','Sepsis2','Sepsis3','Sepsis4')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("./",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                                   #每个基因至少在3个细胞中表达，每一个细胞至少有250个基因表达
                                   min.cells = 3, min.features = 250)
}
#修改名称
names(datalist)=dir_name

#2、细胞质控####
# 批量计算线粒体和rRNA占比
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比
  datalist[[i]] <- sce
  rm(sce)
}
#质控前的
violin=list()
for (i in 1:length(datalist)){
  violin[[i]] <- VlnPlot(datalist[[i]],
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                         pt.size = 0.1, 
                         ncol = 4)
}
pearplot_befor <- CombinePlots(plots = violin , nrow=length(datalist), legend="none")
pearplot_befor
ggsave(filename = 'QC_before.pdf',plot = pearplot_befor,he=15,wi=15)
#样本合并
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
#统计每一个样本的个数

raw_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)

pearplot_befor1<-VlnPlot(sce,
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                        pt.size = 0.1, 
                        ncol = 4)
pearplot_befor1
ggsave(filename = 'QC_before1.pdf',plot = pearplot_befor1,he=7,wi=15)
rm(sce)

#过滤
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset = nFeature_RNA > 200 & 
              nFeature_RNA < 5000 & 
              quantile(percent.mt, 0.98) > percent.mt & percent.mt < 10 &
              quantile(percent.Ribo, 0.99) > percent.Ribo & percent.Ribo > quantile(percent.Ribo, 0.01) & 
              nCount_RNA < quantile(nCount_RNA, 0.97) & nCount_RNA > 1000 )
})
#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)

#过滤前后样本细胞数据的统计
summary_cells <- as.data.frame(cbind(raw_count,clean_count))
counts <- rbind(as.data.frame(cbind(summary_cells[,1],rep("raw",each = length(summary_cells[,1])))),
                as.data.frame(cbind(summary_cells[,2],rep("clean",each = length(summary_cells[,2])))))
counts$sample <- rep(rownames(summary_cells),times =2)
colnames(counts)<- c("count","Stat","sample")
counts[,1] <- as.numeric(counts[,1])
counts$Stat <- factor(counts$Stat, levels=c("raw", "clean"), ordered=TRUE)
fit_cell_count <- ggplot(data =counts, mapping = aes(x = sample, y=count))+ 
  geom_bar(aes(fill = Stat),stat = 'identity', position = 'dodge') + scale_fill_brewer(palette = "Set1") +
  theme(text=element_text(size=10),legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",size = 0.2),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))

fit_cell_count
ggsave(filename = 'fit_cell_count.pdf',plot = fit_cell_count,width = 9,height = 9)

#质控后的小提琴图
violin_after=list()
for (i in 1:length(datalist)){
  violin_after[[i]] <- VlnPlot(datalist[[i]],
                               features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"), 
                               pt.size = 0.1,
                               ncol = 4)
}
pearplot_after <- CombinePlots(plots = violin_after , nrow=length(datalist), legend="none")
pearplot_after
ggsave(filename = 'QC_after.pdf',plot = pearplot_after,he=15,wi=15)
pearplot_after1 <- VlnPlot(sce,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"), 
                         pt.size = 0.1,
                         ncol = 4)
pearplot_after1
ggsave(filename = 'QC_after1.pdf',plot = pearplot_after1,he=7,wi=15)
#质控前后图片的合并
pearplot_befor1
pearplot_after1
qc_merge<- CombinePlots(plots = list(pearplot_befor1,pearplot_after1) , 
                        nrow=2, legend='none')
qc_merge
ggsave(filename = 'qc_merge.pdf',plot = qc_merge,he=9,wi=15)


#保存datalist文件
save(datalist,file = 'datalist.RData')

#3、数据预处理####
#合并数据
library(harmony)
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
rm(datalist)
sce <- NormalizeData(sce, normalization.method = "LogNormalize") 
sce <- FindVariableFeatures(sce)
# #ScaleData矫正什么因素就在regress里面些什么
sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo",
                                                   "percent.mt", "nCount_RNA"), verbose = T)
sce <- RunPCA(sce, verbose = T,npcs = 50)
# 
# install.packages("harmony")
# library(harmony)
sce <- RunHarmony(sce, "orig.ident")
names(sce@reductions)
sce <- RunUMAP(sce,  dims = 1:25, reduction = "harmony")
DimPlot(sce,reduction = "umap",label=T,split.by  = 'orig.ident') 
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:25) 
#harmony_embeddings <- Embeddings(sce.mergeTEN, 'harmony')
#dim(harmony_embeddings)
#p3 <- DimPlot(object = sce.mergeTEN, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
#p4 <- VlnPlot(object = sce.mergeTEN, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
#CombinePlots(plots=list(p3,p4))

#p1+p3

#sce.mergeTEN <- sce.mergeTEN %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  RunTSNE(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)                   
### 可视化前20个高变基因
top20 <- head(VariableFeatures(sce), 20)
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE, size=3.0)

feat_20 <- CombinePlots(plots = list(plot1, plot2),legend="bottom")
feat_20
ggsave(filename = 'feat_20.pdf',plot = feat_20,he=10,wi=15)
#ScaleData
#scale.genes <-  rownames(sce)
#sce <- ScaleData(sce, features = scale.genes)
#样本的分组
meta1<-data.frame(matrix(nrow=length(sce@meta.data$orig.ident), ncol=2)) 
colnames(meta1)=c('Sample','Group1')
meta1$Sample=sce@meta.data$orig.ident
unique(meta1$Sample)
### Group1 Tumor 为原发性肿瘤；Normal：正常
meta1[grep("Normal1",meta1$Sample),]$Group1="Normal"
meta1[grep("Normal2",meta1$Sample),]$Group1="Normal"
meta1[grep("Sepsis1",meta1$Sample),]$Group1="Sepsis"
meta1[grep("Sepsis2",meta1$Sample),]$Group1="Sepsis"
meta1[grep("Sepsis3",meta1$Sample),]$Group1="Sepsis"
meta1[grep("Sepsis4",meta1$Sample),]$Group1="Sepsis"

sce <- AddMetaData(sce, meta1$Sample,col.name = "Sample")
sce <- AddMetaData(sce, meta1$Group1,col.name = "Group1")
save(sce,file = 'sce.RData')



############################################
Dims <- 40
Resolution <- 0.1
sce <- FindNeighbors(object = sce, dims = 1:Dims)
sce <- FindClusters(object = sce, resolution = Resolution)
#颜色
allcolour=c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
              '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69','#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')
length(table(sce@active.ident))
mycolor = allcolour[1:length(table(sce@active.ident))]
#### 按cluster进行占比统计
cluster.frequency.table <- sce@meta.data %>%
  dplyr::count(seurat_clusters) %>%
  dplyr::mutate(freq = n / sum(n)*100) %>%
  ungroup()%>%as.data.frame()

cluster.frequency.table

pie(cluster.frequency.table$n, labels=round(cluster.frequency.table$freq,2),radius=1.0, main = "Percentage of Cluster", col=mycolor)   
legend("right",legend=unique(cluster.frequency.table$seurat_clusters),bty="n",fill=mycolor)

#统计每一个分组，每一个亚群所占的比例
cluster.frequency.sample=data.frame()
for (i in as.character(unique(sce@meta.data$Group1))){
  data1<-sce@meta.data[which(sce@meta.data$Group1==i),]
  dat1 <- data1 %>%
    dplyr::group_by(Group1) %>%
    dplyr::count(seurat_clusters) %>%
    dplyr::mutate(freq = n / sum(n)*100) %>%
    ungroup()%>%as.data.frame()
  cluster.frequency.sample=rbind(cluster.frequency.sample,dat1)
}
head(cluster.frequency.sample)
cluster.freq.sample<-tidyr::spread(data=cluster.frequency.sample[,c("Group1","seurat_clusters","freq")],
                                   key=Group1, value=freq)
cluster.freq.sample[is.na(cluster.freq.sample)]<-0
head(cluster.freq.sample)

#rownames(cluster.freq.sample)=cluster.freq.sample$seurat_clusters
#从内圈到外圈依次是 Normal 、Turmal
cluster.freq<-ggplot(data=cluster.frequency.sample, mapping=aes(x=Group1,y=freq,fill=seurat_clusters))+
  geom_bar(stat='identity',width=0.9)+coord_polar(theta="y",start = 0)+
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank())+
  scale_fill_manual(values=mycolor)
cluster.freq
pdf('cluster_freq.pdf',he=7,wi=9)
cluster.freq
dev.off()
write.csv(cluster.frequency.sample,file ="cluster.frequency.csv")

### UMAP
sce <- RunUMAP(sce, dims=1:Dims, reduction="pca")
###tsne 降维
sce <- RunTSNE(sce, 
               dims=1:Dims, 
               reduction="pca",
               perplexity=30,
               max_iter=1000)

#可视化
sc_umap = DimPlot(sce,cols=mycolor,
                  reduction="umap",
                  #reduction="tsne",
                  label = "T", 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 
sc_umap
ggsave('sc_umap_cluster.pdf',sc_umap,he=7,wi=7)

sc_tsne = DimPlot(sce,cols=mycolor,
                  #reduction="umap",
                  reduction="tsne",
                  label = "T", 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 
sc_tsne
ggsave('sc_tsne_cluster.pdf',sc_tsne,he=7,wi=7)

sc_umap_group1 = DimPlot(sce,cols=mycolor,group.by='Sample',
                         reduction="tsne",
                         label = "T", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 

sc_umap_group1
ggsave('sc_tsne_sample.pdf',sc_umap_group1,he=7,wi=7)

sc_umap_group2 = DimPlot(sce,cols=mycolor,group.by='Group1',
                         reduction="tsne",
                         label = "T", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 

sc_umap_group2
ggsave('sc_tsne_group.pdf',sc_umap_group2,he=7,wi=7)

###############################################
Logfc = 0.25
#差异基因时最小的表达比例
Minpct = 0.25
DefaultAssay(sce) <- "RNA"
sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, 
                              min.pct = Minpct,only.pos = T)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
length(unique(sce.markers$gene))
head(sce.markers)

write.table(sce.markers,'scRNA_marker_gene.txt',quote = F,row.names = F,sep='\t')
### 选择前5个marker基因
Top5 <- sce.markers %>% 
  group_by(cluster) %>% 
  slice_max(n =5, order_by = avg_log2FC)  
Top5 <- unique(Top5$gene)

sc_marker_dotplot <- DotPlot(object = sce, 
                             features = Top5,
                             cols=c("blue", "red"),
                             scale = T)+ 
  RotatedAxis()+ ggtitle("Top 5 Marker Genes")+ 
  theme(plot.title = element_text(hjust = 0.5)) 

sc_marker_dotplot
ggsave(filename = 'sc_marker_dotplot.pdf',
       plot = sc_marker_dotplot,
       height = 9,width = 25)
#热图展示
library(viridisLite)
sc_marker_heatmap<- DoHeatmap(object = sce,
                              features = Top5,
                              group.colors = mycolor,
                              label = F) + 
  ggtitle("Top 5 Marker Genes") + 
  theme(plot.title = element_text(hjust = 0.5)) 
sc_marker_heatmap
ggsave(filename = 'sc_marker_heatmap.pdf',
       plot = sc_marker_heatmap,
       width = 12,height = 12)

save(sce,file = 'sce2.RData')


#############################################################
markers <-  c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',  ## T 
                   'CD19', 'CD79A', 'MS4A1' ,  ## B 
                   'IGHG1', 'MZB1', 'SDC1', ## Plasma 
                   'CD68', 'CD163', 'CD14', 
                   'VCAN', 'FCN1', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'C1QA',  'C1QB',  # mac
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'MKI67' , 'TOP2A','LYZ', ##myloid
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'CLEC10A', 'CLEC12A','LILRA4', ##pDC
                   'FGFBP2','NCR1','KLRK1','KLRC1','NCAM1','FCGR3A', # NK 
                   'KIT', ##MAST
                   'CSF3R','ITGAM','ITGB2','FUT4','FCGR3B','FCGR2A', ##Neutrophil
                   'FGF7','MME', 'ACTA2','VWF','COL1A1', ## human  fibo 
                   'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo 
                   'PECAM1', 'VWF',  ## endo 
                   'EPCAM' , 'KRT19','KRT7' # epi 
)
alldata <- ScaleData(sce, 
                     features = markers, 
                     assay = "RNA")

DoHeatmap(alldata, 
          features = markers,
          group.by = "seurat_clusters",
          assay = "RNA")


library(SingleR)
library(celldex)
library(BiocParallel)
library(Seurat)

###下载
hpca.se <- HumanPrimaryCellAtlasData()
#load('sce2.RData')
#获取基因的表达谱的count数据
testdata <- GetAssayData(sce, slot="data")
#获取聚类的亚群
clusters <- sce@meta.data$seurat_clusters
pred.sce <- SingleR(test =  testdata, 
                    ref = hpca.se, 
                    labels = hpca.se$label.fine,
                    method = "clusters",
                    clusters = clusters, 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

plotScoreHeatmap(pred.sce)

celltype = data.frame(ClusterID=rownames(pred.sce), celltype=pred.sce$labels
                      , stringsAsFactors = F)
celltype
write.table(celltype,'celltype.txt',quote = F,sep = '\t',row.names = F)

#修改亚群的名称
sce <- RenameIdents(object = sce, 
                      "0"="T_cell", 
                     "1"="Neutrophil", 
                     "2"="T_cell", 
                     "3"="B_cell", 
                     "4"="Macrophage", 
                     "5"="Epithelial", 
                     "6"="Plasma", 
                     "7"="Fibroblasts", 
                     "8"="T_cell", 
                     "9"="Endothelial", 
                     "10"="Mast_cell")
                     

length(table(sce@active.ident))    
sce$celltype <- sce@active.ident
DimPlot(sce, group.by = "celltype")
   allcolour=c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
              '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69','#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')
#length(table(sce@active.ident))
mycolor = allcolour[1:length(table(sce@active.ident))]

umap_celltype<-DimPlot(sce,cols=mycolor,
                       reduction="tsne",
                       label = "T", 
                       pt.size = 0.5,
                       label.size = 5)

umap_celltype
ggplot2::ggsave('tsne_celltype.pdf',plot = umap_celltype,he=7,wi=9)
save(sce, file = 'sce3.RData')

#从这里开始---------------------------------
rm(list=ls())
load('../1.单细胞/Necroptosis.sce.RData')
sce
head(sce@meta.data)
table(sce$Group1, sce$celltype)

markers <- c('PTPRC','CD3D','CD3E','CD4','CD8A','CD19', 'CD79A', 'MS4A1' ,'IGHG1', 'MZB1', 'SDC1','TPSAB1','TPSB2',
'KIT','RCVRN','FPR1','ITGAM','C1QA', 'C1QB','CSF3R','ITGB2','FCGR3B','FCGR2A','FGF7','MME', 'ACTA2','VWF','COL1A1','PECAM1', 'EPCAM' ,'KRT19','KRT7')
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
              '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69','#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')

sce$celltype1 <- paste(sce$Necroptosis.score.group, sce$celltype, sep = "_") 
DotPlot(sce, features = c("SOCS3","PIM3","MARCKS","CEBPD","CEBPB","JUNB"), 
        group.by = "celltype"  ) #+ coord_flip()

table(sce$celltype1)

saveRDS(sce, "GSE167363_scRNA.rds")
#install.packages("remotes")  
#remotes::install_github("lyc-1995/MySeuratWrappers")#通过链接安装包  
library(MySeuratWrappers)
library(ggplot2)

DefaultAssay(sce) <- "RNA"
Idents(sce) <- "celltype"
VlnPlot(sce, features = markers,
        stacked=T,pt.size=0,
        cols = my36colors,#颜色
        direction = "horizontal", #水平作图
        x.lab = '', y.lab = '')+#横纵轴不标记任何东西
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())#不显示坐标刻度  
        
VlnPlot(sce, features = markers, 
        group.by = 'celltype',stack = T,cols = my36colors,#颜色
        pt.size = 0)&
  theme(legend.position = "none",
        panel.spacing = unit(x = 0, units = 'lines'),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.ticks.length = unit(x = 0, units = 'cm'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
ggplot2::ggsave('VlnPlot.pdf',he=12,wi=8)       

#气泡热图3整合修饰------------------------------------------------
#install.packages("devtools")
#devtools::install_github("Simon-Leonard/FlexDotPlot")
library(FlexDotPlot)
library(Seurat)
library(ggplot2)
dp = DotPlot(sce, features = markers) + RotatedAxis()
dot_plot(dp$data[,c(3,4,1,2,5)],
         size_var = "pct.exp", 
         col_var = "avg.exp.scaled",
         size_legend = "Percent Expressed", 
         col_legend = "Average Expression",
         x.lab.pos = "bottom", 
         y.lab.pos = "right",
         display_max_sizes = F, 
         shape.scale=6,
         hclust_method = "ward.D2",
         dend_x_var = c("pct.exp", "avg.exp.scaled"),
         dend_y_var = c("pct.exp", "avg.exp.scaled"),
         text.size = 0.5,
         text.vjust = 0.5,
         size.breaks.number=6,
         color.breaks.number=4,
         x.lab.size.factor = 0.8,
         cols.use = c("#1e4668", "#59829e","white","#82677e" ,"#492952" ),
         y.lab.size.factor = 1)
#ggsave('气泡.pdf',he=8,wi=12)          
        
#Necroptosis-score计算免疫原性死亡相关评分----------------------------------
Necroptosis.gene=read.delim('Necroptosis.txt',sep='\t',header = F)
Necroptosis.gene=Necroptosis.gene$V1
Necroptosis.gene
Necroptosis.gene=list(Necroptosis.gene)
names(Necroptosis.gene)='Necroptosis'
ssGSEA_score<-function(dat, kcdf = c("Gaussian", "Poisson")[1], GeneSet =  NULL) {
  # dat：行为基因，列为样本
  # kcdf：log2(FPKM、TPM)数据使用 Gaussian；Counts 数据使用 Poisson
  library(GSVA)

  if (is.null(GeneSet)) {
    gene_set <- read.table(paste0(basedir,"/SourceFiles/PMID_28052254.txt"), header = T, stringsAsFactors = F, sep = '\t')[, 1:2]
    gene_set <- split(as.matrix(gene_set)[,1], gene_set[,2])
  } else {
    gene_set <- GeneSet
  }

  ssgsea_res <- gsva(as.matrix(dat),
                     gene_set,
                     method='ssgsea',
                     kcdf=kcdf,
                     abs.ranking=TRUE)
  ssgsea_res <- as.data.frame(t(ssgsea_res))
  return(ssgsea_res)
}

sce.Necroptosis_score <- ssGSEA_score(dat = as.matrix(sce@assays$RNA@counts),kcdf='Poisson',
                              GeneSet = Necroptosis.gene)
head(sce.Necroptosis_score)
sce.Necroptosis_score=sce.Necroptosis_score[rownames(sce@meta.data),,drop=F]
sce$Necroptosis.score=sce.Necroptosis_score[,1]

#TSNE的分布
color <- c('lightgrey', 'blue','seagreen2')#设置颜色 
fig1c<-FeaturePlot(object = sce, features = "Necroptosis.score",reduction = 'tsne',cols = color, pt.size = 1)
fig1c
ggsave('fig1c.pdf',fig1c,height = 6,width = 7)

#各个亚群的分布
fig1a<-DimPlot(sce,group.by = 'celltype',cols=my36colors,
                 reduction="tsne", 
                 label = "T", pt.size = 0.5,
                 label.size = 5) 
fig1a
ggsave('fig1a.pdf',fig1a,height = 6,width = 9)

#各个亚群表达的Necroptosis的小提琴图
Necroptosis.score.cell=sce@meta.data[,c("celltype","Necroptosis.score")]
head(Necroptosis.score.cell)
library(ggplot2)
fig1d <- ggplot(Necroptosis.score.cell,aes(x=celltype,y=Necroptosis.score))+
         geom_violin(aes(fill=celltype))+theme_bw()+scale_fill_manual(values = my36colors)
         theme(legend.position = 'none',
         axis.text.x = element_text(angle = 30,hjust = 1))+
         xlab('')
fig1d
ggsave('fig1d.pdf',fig1d,height = 7,width = 9)

#高低Necroptosis-score细胞的差异分析
sce$Necroptosis.score.group=ifelse(sce$Necroptosis.score>median(sce$Necroptosis.score),'High','Low')
DefaultAssay(sce) <- "RNA"
Idents(sce)<-'Necroptosis.score.group'
Necroptosis.markers <- FindAllMarkers(object = sce,logfc.threshold = 0, min.pct = 0.35,only.pos = T)
head(Necroptosis.markers)
table(Necroptosis.markers$cluster)
write.table(Necroptosis.markers,'Necroptosis.markers.txt',quote = F,row.names = F,sep='\t')
save(sce,file = 'Necroptosis.sce.RData') 

library(Seurat)
library(dplyr)
library(reshape2)
library(plyr)
library(ggplot2)
source("./Singlecellratio_plotstat.R")

sce@meta.data$Group = sce@meta.data$Group1
#sce$celltype<-Idents(sce)
my_comparisons <- list(c("Normal", "Sepsis"))

#先做一个柱状图
pdf("比例.pdf",width = 6,height = 8)
Singlecellratio_plotstat(sce, group_by = "Group",
                         meta.include = c("Group","orig.ident"),
                         color_by = 'celltype')
dev.off()

#分组箱线图1
Singlecellratio_plotstat(sce, group_by = "Group",
                         meta.include = c("Group","orig.ident"),
                         comparisons = my_comparisons, color_by = 'Group',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =3)

#分组箱线图2
Singlecellratio_plotstat(sce, group_by = "Group",
                         meta.include = c("Group","orig.ident"),
                         comparisons = my_comparisons, color_by = 'orig.ident',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =3)

#分组点不同形状
pdf("分组.pdf",width = 12,height = 9)
Singlecellratio_plotstat(sce, group_by = "Group",
                         meta.include = c("Group","orig.ident"),
                         comparisons = my_comparisons, color_by = 'orig.ident',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =4,shape_by = 'Group')
dev.off()                         


pbmc <- readRDS("ICC_seu_wenzai.rds")
head(pbmc@meta.data)
table(pbmc$Sample)
DimPlot(pbmc, group.by = "Sample")
median(pbmc$nFeature_RNA)
