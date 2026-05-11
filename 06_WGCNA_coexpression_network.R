
dir.create('results')
library(WGCNA)
tcga.dat<-read.delim('../2.数据收集/datExpr.txt',sep='\t',header = T,check.names = F,row.names = 1)
Necroptosis.gene=read.delim('../1.单细胞/Necroptosis.txt',sep='\t',header = F)
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

tcga.Necroptosis <- ssGSEA_score(dat = as.matrix(tcga.dat),kcdf='Gaussian',
                              GeneSet = Necroptosis.gene)
head(tcga.Necroptosis)
Necroptosis.related.gene=read.delim('../1.单细胞/Necroptosis.markers.txt',sep='\t',header = T)
table(Necroptosis.related.gene$cluster)
#
gse_imcell_wgcna <- na.omit(tcga.dat[unique(Necroptosis.related.gene$gene),])
gse_imcell_wgcna[1:4,1:4]
#选择有意义的亚群（从03这一步来的，选择肿瘤发生，预后等相关的）进行后期的WGCNA分析，计算该亚群与module之间的相关性
gse_imcell_cli <- tcga.Necroptosis

com.sample=intersect(rownames(gse_imcell_cli),colnames(gse_imcell_wgcna))
range(gse_imcell_wgcna)
gse_imcell_wgcna=gse_imcell_wgcna[apply(gse_imcell_wgcna,1,function(x){return(sd(x)>0)}),]


datExpr=t(gse_imcell_wgcna)
gse_imcell_cli=gse_imcell_cli[com.sample,,drop=F]
range(datExpr)
head(gse_imcell_cli)
gsg = goodSamplesGenes(datExpr, verbose = 3)
#这一步一定要是TURE，如果是FALSE要重新过滤数据
gsg$allOK


dim(datExpr)

sampleTree = hclust(dist(datExpr), method = "complete")
plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")
# abline(h=140,col='red')

fpkmSamples = rownames(datExpr) 
traitSamples =rownames(gse_imcell_cli) 
traitRows = match(fpkmSamples, traitSamples) 
datTraits = gse_imcell_cli[traitRows,,drop=F] 
rownames(datTraits) 
collectGarbage()
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度

pdf(file="results/Fig2a.pdf",width=20,height=12) 
plotDendroAndColors(sampleTree, traitColors, 
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap",cex.colorLabels = 1.5, cex.dendroLabels = 1, cex.rowText = 2) 
dev.off() 


clust = cutreeStatic(sampleTree, cutHeight = 60000, minSize = 3)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, 
                        powerVector = powers, 
                        verbose = 3)


cex1 = 0.85
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()
softPower=sft$powerEstimate
softPower
#5
net = blockwiseModules(datExpr, power = softPower, maxBlockSize = 5000,
                       TOMType = "unsigned", minModuleSize = 20,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "TPM-TOM-40",
                       deepSplit = 3,
                       verbose = 3)
table(net$colors)
mergedColors = labels2colors(net$colors)
table(mergedColors)
length(table(mergedColors))
pdf('results/fig2b.pdf',width = 8,height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    groupLabels = c("Module colors","GS.weight"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


MEs_col = net$MEs[rownames(gse_imcell_cli), ]
colnames(MEs_col) = labels2colors(as.numeric(gsub('ME','',colnames(net$MEs))))

spms=gse_imcell_cli[, unique(colnames(gse_imcell_cli)),drop=F]

dim(spms)[1]
modTraitCor = cor(MEs_col, as.data.frame(spms), use = "p")
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf('results/Fig2c.pdf',width = 5,height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(spms), 
               yLabels = colnames(MEs_col), 
               cex.lab = 1, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

MEs = moduleEigengenes(datExpr[rownames(gse_imcell_cli),], mergedColors)$eigengenes
MEs = orderMEs(MEs)
colnames(MEs)=gsub('^ME','',colnames(MEs))
me.inds=which(colnames(MEs)!='grey')
if(length(me.inds)==1){
  cns=colnames(MEs)[me.inds]
  MEs=as.data.frame(MEs[,me.inds])
  colnames(MEs)=cns
}else if(length(me.inds)>1){
  MEs=MEs[,me.inds]
}

nGenes = ncol(datExpr)
nSamples = nrow(datExpr[rownames(gse_imcell_cli),])
geneTraitSignificance = as.data.frame(cor(datExpr[rownames(gse_imcell_cli),],spms, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(MEs)));
names(geneTraitSignificance) = paste("GS.", colnames(spms), sep="");
names(GSPvalue) = paste("p.GS.", colnames(spms), sep="")
geneModuleMembership = as.data.frame(cor(datExpr[rownames(gse_imcell_cli),], MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", names(geneModuleMembership), sep="");
names(MMPvalue) = paste("p.MM", names(geneModuleMembership), sep="");
write.table(cbind(Genes=row.names(geneTraitSignificance),geneTraitSignificance,GSPvalue,Module=mergedColors)
            ,file = "results/GeneTraitSignificance.txt",sep = '\t',quote = F,row.names = F)

#salmon 模块####
table(mergedColors)
salmon_genes <- colnames(datExpr)[which(mergedColors=='salmon')]
#富集分析
length(salmon_genes)#140

write.table(salmon_genes,'results/module.gene.txt',quote = F,row.names = F,sep='\t',col.names = F)
#
pdf('results/fig2d.pdf',he=7,wi=7)
verboseScatterplot(abs(geneModuleMembership[salmon_genes,"MMsalmon"]),
                         abs(geneTraitSignificance[salmon_genes,"GS.Necroptosis"]),
                         xlab = paste0("Module Membership in salmon module"),
                         ylab = "Gene significance for Necroptosis score",
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'salmon')
dev.off()  
