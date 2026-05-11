
dir.create('results')

diffSigOut=rbind(id=colnames(diffSig),diffSig)
tcga.diff.fit=diffSigOut
module.gene=read.delim('../3.WGCNA/results/module.gene.txt',sep='\t',header = F)
module.gene=module.gene$V1
module.gene
com.gene=intersect(module.gene,rownames(tcga.diff.fit))
length(com.gene)
#venn图
venn_plot=function(data_list,fill=NULL,ven=T){
  if(length(data_list)<=4&ven){
    library(VennDiagram)
    if(is.null(fill)){
      fill=c('red', 'blue','green','yellow')[1:length(data_list)]
    }
    g=venn.diagram(data_list, filename = NULL,margin = 0.2, 
                   fill = fill, alpha = 0.50, col = 'black', cex = 1, fontfamily = 'serif'
                   #,cat.col = c('black', 'black', 'black', 'black')
                   ,cat.cex = 1, cat.fontfamily = 'serif')
    grid.draw(g)
  }else if(length(data_list)<=7&ven){
    library(venn)
    if(is.null(fill)){
      fill=mg_colors[1:length(data_list)]
    }
    g=venn(data_list, zcolor = fill,ggplot=F,box=F) 
  }else{
    anm=unique(unlist(data_list))
    tbs=cbind()
    for(i in 1:length(data_list)){
      a1=rep(0,length(anm))
      a1[anm%in%data_list[[i]]]=1
      tbs=cbind(tbs,a1)
    }
    colnames(tbs)=names(data_list)
    row.names(tbs)=anm
    tbs=crbind2DataFrame(tbs)
    g=UpSetR::upset(tbs, nsets = ncol(tbs), nintersects = 30, mb.ratio = c(0.5, 0.5),
                    order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE)
    )
    return(g)
  }
}
pdf('results/fig2f.pdf',height = 5,width = 6)
venn_plot(data_list = list(module.gene=module.gene,
                           deg=rownames(tcga.diff.fit)),fill=c('#E5D2DD', '#53A85F'),ven = T)
dev.off()
write.table(com.gene,'results/com.gene.txt',quote = F,row.names = F,sep='\t')
