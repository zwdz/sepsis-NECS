
dir.create('results')
tcga.dat=read.table('../2.GEO数据收集/treatData.txt',sep='\t',header = T,check.names = F,row.names = 1)
tcga.dat[1:5,1:5];dim(tcga.dat)
tcga.risk.cli<-read.table('../7.immune/Modelrisk.txt',sep='\t',header = T,check.names = F)


#免疫相关通路
immnu_patwhays<-list.files('../7.immune/immnue_pathway/',pattern = '.gmt$')
immnu_patwhay.gene= as.data.frame(data.table::rbindlist(lapply(immnu_patwhays, function(x){
  print(x)
  a=clusterProfiler::read.gmt(paste0('../7.immune/immnue_pathway/',x))
  return(a)
})))
head(immnu_patwhay.gene)
immnu_patwhay.gene$term=gsub('KEGG_','',immnu_patwhay.gene$term)
immnu_patwhay.gene=split(x = immnu_patwhay.gene$gene,f = immnu_patwhay.gene$term)
immnu_patwhay.gene
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

imm.riskscore <- ssGSEA_score(dat = tcga.dat,kcdf='Gaussian',
                              GeneSet = immnu_patwhay.gene)
head(imm.riskscore)
write.table(imm.riskscore,'results/imm.riskscore.txt',quote = F,row.names = T,sep='\t')

#计算p值
diff_pathway<-function(dat,group){
  dat=data.frame(cluster=group,t(dat))
  gr=as.character(unique(group))
  dat1=t(dat[dat$cluster==gr[1],-1])
  dat2=t(dat[dat$cluster==gr[2],-1])
  pathway=unique(c(rownames(dat1),rownames(dat2)))
  p_vale=data.frame()
  for (i in pathway){
    dd1=t.test(as.numeric(dat1[i,]),as.numeric(dat2[i,]))$p.value
    p_vale=rbind(p_vale,data.frame(pathway=i,p.value=dd1))
  }
  return(p_vale)
}
pathwy_p<-diff_pathway(dat=t(imm.riskscore[tcga.risk.cli$ID,]),group=tcga.risk.cli$risk)
head(pathwy_p)
pathwy_p=pathwy_p[pathwy_p$p.value<0.05,]
pathwy_p$lab=paste0(pathwy_p$pathway,ifelse(pathwy_p$p.value<0.001,'***',ifelse(pathwy_p$p.value<0.01,'**',ifelse(pathwy_p$p.value<0.05,'*',''))))
pathwy_p$lab
imm.riskscore[1:4,1:4]
imm.riskscore1=merge(pathwy_p[,c("pathway","lab")],t(imm.riskscore),by.y=0,by.x='pathway')
head(imm.riskscore1)
rownames(imm.riskscore1)=imm.riskscore1$lab
imm.riskscore1=imm.riskscore1[,-c(1,2)]
imm.riskscore1[1:4,1:5]

library(pheatmap)
anno_col=data.frame(risk=tcga.risk.cli$risk,row.names = tcga.risk.cli$ID)
anno_col=anno_col[order(anno_col$risk),,drop=F]
pdf('fig9b1.pdf',height = 7,width = 15)
pheatmap(mat = imm.riskscore1[,rownames(anno_col)],scale = 'row',
         annotation_col = anno_col,
         cluster_cols = F, cluster_rows = T,
         show_rownames = T, show_colnames = F,
        color = colorRampPalette(c("blue", "#FFFFFF", "red"))(100),
         annotation_colors = list(risk = c('high' = "red", 'low' = "blue")))
dev.off()


#CIBERSORT
tcga.cibersort <- deconvo_tme(eset = tcga.dat,
                              method = "cibersort",
                              arrays = FALSE,
                              perm = 100) #
save(tcga.cibersort,file = '../7.immune/results/tcga.cibersort.RData')

load('../7.immune/results/tcga.cibersort.RData')
tcga.cibersort=as.data.frame(tcga.cibersort)
tcga.cibersort[1:4,1:4]
tcga.cibersort=tcga.cibersort[,1:23]
vioplot_plot=function(rt,normal,tumor,normal.name,tumor.name){
  library(vioplot)
  rt=rt[c(normal,tumor),]
  normal=length(normal)
  tumor=length(tumor)
  par(las=1,mar=c(10,6,3,3))
  x=c(1:ncol(rt))
  y=c(1:ncol(rt))
  plot(x,y,
       xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
       main="",xlab="", ylab="Fraction",
       pch=21,
       col="white",
       xaxt="n")
  for(i in 1:ncol(rt)){
    if(sd(rt[1:normal,i])==0){
      rt[1,i]=0.001
    }
    if(sd(rt[(normal+1):(normal+tumor),i])==0){
      rt[(normal+1),i]=0.001
    }
    normalData=rt[1:normal,i]
    tumorData=rt[(normal+1):(normal+tumor),i]
    vioplot(normalData,at=3*(i-1),lty=1,add = T,col = "#1e4668")
    vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = "#990033")
    wilcoxTest=wilcox.test(normalData,tumorData)
    p=round(wilcoxTest$p.value,3)
    mx=max(c(normalData,tumorData))
    lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
    #text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.05, paste0("p=",p),''), cex = 0.8)
    text(x=3*(i-1)+0.5, y=mx+0.02, labels= paste0("p=",p), cex = 0.8)
  }
  text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
  legend(x = "topleft", box.col = "black",
         bg ="white", box.lwd = 0.5 , cex = 0.5,
         legend=c(normal.name, tumor.name), 
         fill = c("#1e4668","#990033"))
}
high.sample=tcga.risk.cli[which(tcga.risk.cli$risk=='high'),"ID"]
low.sample=tcga.risk.cli[which(tcga.risk.cli$risk=='low'),"ID"]
rownames(tcga.cibersort)=tcga.cibersort$ID
colnames(tcga.cibersort)=gsub('_CIBERSORT','',colnames(tcga.cibersort))
pdf('../7.immune/results/fig9c.pdf',height = 6,width = 10)
vioplot_plot(rt = tcga.cibersort[,-1],normal = low.sample,tumor =high.sample,normal.name = 'Risk-low',tumor.name = 'Risk-high')
dev.off()

#22种免疫细胞评分与关键基因的相关性分析
sig.gene=c("CEBPD","CEBPB","MARCKS","SOCS3","PIM3","JUNB")
sig.gene
tcga.gene.exp=t(tcga.dat[sig.gene,])
tcga.gene.exp[1:4,1:4]
ccor=Hmisc::rcorr(x=as.matrix(tcga.gene.exp),
                  y=as.matrix(tcga.cibersort[,-1]),
                  type = 'spearman')

ccor.r=ccor$r[sig.gene,colnames(tcga.cibersort)[-1]]
ccor.p=ccor$P[sig.gene,colnames(tcga.cibersort)[-1]]



pdf('results/fig9e.pdf',height = 9,width = 8)
corrplot::corrplot(as.matrix(t(ccor.r)), 
                   p.mat = as.matrix(t(ccor.p)),
                   mar = c(0,0,1,1),
                   col=colorRampPalette(c("#003366","white","#990033"))(100),
                   tl.srt = 90,
                   tl.cex = 1,
                   tl.col = 'black',
                   tl.offset = 0.5,
                   cl.pos = c("b","r","n")[1], 
                   cl.align.text = 'l',
                   cl.length = 5,
                   cl.ratio = 0.1,
                   cl.cex = 0.8,
                   addgrid.col = 'white',
                   method = 'color',
                   insig = 'label_sig',
                   sig.level=c(0.001,0.01,0.05),
                   pch.cex=1,
                   is.corr=T,
                   xpd=T
)
dev.off()

##22种免疫细胞与风险得分的相关性分析
#tcga.risk.cli<-read.delim('../09.KM/results/tcga.risk.cli.txt',sep='\t',header = T)
rownames(tcga.risk.cli)=tcga.risk.cli$ID
tcga.risk.cli=tcga.risk.cli[tcga.cibersort$ID,]
imm.risk.ccor=Hmisc::rcorr(as.matrix(cbind.data.frame(tcga.cibersort[,-1],
                                     Riskscore=tcga.risk.cli$riskscore)),
                  type = 'spearman')
imm.risk.ccor.r=reshape2::melt(imm.risk.ccor$r)
imm.risk.ccor.p=reshape2::melt(imm.risk.ccor$P)
imm.risk.ccor.r=imm.risk.ccor.r[which(imm.risk.ccor.r$Var1 =='Riskscore' & imm.risk.ccor.r$Var2 !='Riskscore'),]
imm.risk.ccor.p=imm.risk.ccor.p[which(imm.risk.ccor.p$Var1 =='Riskscore' & imm.risk.ccor.p$Var2 !='Riskscore'),]

imm.risk.ccor.res=merge(imm.risk.ccor.r,imm.risk.ccor.p,by=c('Var1','Var2'))
head(imm.risk.ccor.res)
colnames(imm.risk.ccor.res)=c('Riskscore','cell','cor','p')
fig9f<-ggplot(imm.risk.ccor.res,aes(cor,forcats::fct_reorder(cell,cor))) +
  geom_segment(aes(xend=0,yend=cell)) +
  geom_point(aes(col=p,size=abs(cor))) +
  scale_colour_gradientn(colours=c("#F9791EFF","#017A4AFF")) +
  scale_size_continuous(range =c(2,8))  +
  theme_classic() +
  ylab(NULL)
fig9f
ggsave('results/fig9f.pdf',fig9f,height = 7,width = 8)

