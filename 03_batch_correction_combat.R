
library(TCGAbiolinks)
library(sva)
library(cluster)
library(oompaBase)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# GSE65858-----------------------
gse65858.expr <- read.delim("gse65858.expr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gse65858.expr[1:3,1:3] # 检查一下，行名很干净，不用处理
gse65858.expr <- as.data.frame(round(gse65858.expr,2))
write.csv(gse65858.expr, "easy_input_GSE65858.csv", quote = F)


# 检查两类数据的量级
#通过取对数，让两类数据的数量级达到同一范围。
library(limma)               #引用包
expFile="GSE65682.txt"     #表达数据文件
#conFile="s1.txt"             #对照组样品文件
#treatFile="s2.txt"           #实验组样品文件
#setwd("C:\\biowolf\\m6A\\06.normalize")      #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
write.table(data, file="normalize.txt", sep="\t", quote=F, col.names=T)

# TCGA-----------------
tcga.expr <- read.table("normalize.txt", row.names = 1,sep="\t",header=T,check.names=F)
#tpms.hugo[1:3, 1:3]
#range(tpms.hugo) 
# 量级很大，所以取对数
#tcga.expr <-tpms.hugo
range(tcga.expr) # 取完对数后，量级在0-20

# GEO microarray-------------------------------
gse54514.expr <-  read.table("GSE54514.txt", row.names = 1,sep="\t",header=T,check.names=F)
gse13904.expr <- read.table("GSE13904.txt", row.names = 1,sep="\t",header=T,check.names=F)
gse95233.expr <- read.table("GSE95233.txt", row.names = 1,sep="\t",header=T,check.names=F)

range(gse54514.expr) # 量级在0-20
range(gse13904.expr)
range(gse95233.expr)
 # 量级在0-20

# 合并数据集、检查批次效应
#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0'
# 设置颜色
red<- '#53A85F'
blue <- "#E5D2DD"
yellow <- '#F1BB72'
green <- '#F3B1A0'

# 提取在三套数据中都出现的基因
comgene <- intersect(intersect(rownames(gse54514.expr), rownames(gse13904.expr)), rownames(gse95233.expr))
comgene <- intersect(comgene, rownames(tcga.expr))
# 合并三套数据
combined.expr <- cbind.data.frame(tcga.expr[comgene,],
                                  gse54514.expr[comgene,],
                                  gse13904.expr[comgene,],
                                  gse95233.expr[comgene,])

# 绘制PCA散点图，检查批次效应
batchPCA(indata = t(scale(t(combined.expr))),
         batch = rep(c("GSE65682","GSE54514","GSE13904","GSE95233"), times = c(ncol(tcga.expr),ncol(gse54514.expr),ncol(gse13904.expr),ncol(gse95233.expr))),
         fig.dir = ".",
         PCA.fig.title = "Raw PCA for combined expression profile",
         cols = c(red,blue, yellow, green),
         showID = F,
         cex = 1,
         showLegend = T) # 可以看到三个数据集(批次)完全分开，说明有很严重的批次效应
range(combined.expr)


#Raw PCA for combined expression profile.pdf

# combat去除批次效应，再次检查
# 去除批次效应
batch <- data.frame(batch = rep(c("GSE65682","GSE54514","GSE13904","GSE95233"), times = c(ncol(tcga.expr),ncol(gse54514.expr),ncol(gse13904.expr),ncol(gse95233.expr))))
modcombat = model.matrix(~1, data = batch)
combined.expr.combat <- as.data.frame(ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))

# 输出到文件
# 这里保存前三个基因，便于传输和查看
write.csv(combined.expr.combat[1:3,], "output_combined_expr.csv", quote = F)
# 实际使用时请运行下面这行，获得全部基因在三个数据集中的表达矩阵
#write.csv(combined.expr.combat, "output_combined_expr.csv", quote = F)

# 绘制PCA散点图，检查批次效应
batchPCA(indata = t(scale(t(combined.expr.combat))),
         batch = rep(c("GSE65682","GSE54514","GSE13904","GSE95233"), times = c(ncol(tcga.expr),ncol(gse54514.expr),ncol(gse13904.expr),ncol(gse95233.expr))),
         fig.dir = ".",
         PCA.fig.title = "Combat PCA for combined expression profile",
         cols = c(red,blue, yellow, green),
         showID = F,
         cex = 1,
         showLegend = T) # 可以看到三个数据集(批次)混杂在了一次，说明批次效应被基本消除
range(combined.expr.combat)
