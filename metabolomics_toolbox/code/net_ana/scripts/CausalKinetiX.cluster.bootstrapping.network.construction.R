# network construction and clustering based on bootstraping samples
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(R.matlab)
require(ggplot2)
require(dplyr)
require(RCy3)
require(ComplexHeatmap)
require(circlize)

# the user will need to modify the path here for local run
comp="/Users/yuewu/"
pardir=paste0(comp,"Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/");
# OR
# pardir=paste0(getwd(),"/")
datadir=paste0(pardir,"result_data/")
bootstrapdir=paste0(pardir,"result/clust_bootstrapping/result/")
resdir=paste0(pardir,"result/clust_bootstrapping/")

strlist=c("_peakwise")
nbootstrp=100
# open cytoscape first before the following
cytoscapePing()
cytoscapeVersionInfo()
# set the default setting in cytoscape of clustering to "Create new clustered network"=TRUE, "Restore inter-cluster edges after layout"=TRUE
listtab=vector(mode="list")
for(add_str in strlist){
  matres=readMat(paste0(datadir,"CausalKinetiX_input",add_str,".mat"))
  compds=unlist(matres[["namesall"]])
  ncompd=length(compds)
  counttab=matrix(0,nrow=ncompd,ncol=ncompd)
  rownames(counttab)=compds
  colnames(counttab)=compds
  for(ibootstrap in 1:nbootstrp){
    filepath=paste0(bootstrapdir,"res",add_str,"/",ibootstrap,"/net.CausalKinetiX",add_str,".txt")
    edgetab=read.table(filepath,sep="\t",header=TRUE)
    createNetworkFromDataFrames(edge=edgetab,title="causalkinetix_network",collection="causalkinetix_network")
    clustermaker_url<-paste0("cluster glay undirectedEdges=true createGroups=false")
    commandsGET(clustermaker_url)
    clusttab=getTableColumns(table="node")
    subtab=clusttab[,c("name","__glayCluster")]
    clusters=unique(subtab[,"__glayCluster"])
    for(cluster in clusters){
      clut_mem=subtab[subtab[,"__glayCluster"]==cluster,"name"]
      ind=which(compds%in%clut_mem)
      counttab[ind,ind]=counttab[ind,ind]+1
    }
    deleteAllNetworks()
  }
  listtab[[paste0("res",add_str)]]=counttab
}
save(listtab,file=paste0(resdir,"connection_tab.RData"))

# Heatmap cluster plot for each peak
set.seed(1)
connecmat=listtab[["res_peakwise"]]
ratiomat=connecmat/nbootstrp
diag(ratiomat)=NA
title="peak clustering frequency"
compdlabels=colnames(ratiomat)
colormap=colorRamp2(c(0,0.7),c("white","red"))
compdlabels=str_replace_all(string=compdlabels,pattern="\\s*[-\\d\\.]+$",replacement="")
ha=HeatmapAnnotation(df=as.data.frame(compdlabels))
pdf(paste0(resdir,"peak_cluster_frequency.pdf"),width=14,height=7)
Heatmap(ratiomat,col=colormap,
        cluster_rows=FALSE,cluster_columns=FALSE,
        show_row_names=FALSE,show_column_names=FALSE,
        name=title,column_title=title,
        show_heatmap_legend=TRUE,
        column_title_gp=gpar(fontsize=12),
        row_names_gp=gpar(fontsize=2),
        top_annotation=ha)
dev.off()

# Heatmap cluster plot for each peak clustered
set.seed(1)
connecmat=listtab[["res_peakwise"]]
ratiomat=connecmat/nbootstrp
diag(ratiomat)=NA
title="peak clustering frequency clustered"
compdlabels=colnames(ratiomat)
compdlabels=str_replace_all(string=compdlabels,pattern="\\s*[-\\d\\.]+$",replacement="")
pdf(paste0(resdir,"peak_cluster_frequency_clustered.pdf"),width=14,height=7)
Heatmap(ratiomat,col=colormap,
        cluster_rows=TRUE,cluster_columns=TRUE,
        show_row_names=TRUE,show_column_names=FALSE,
        name=title,column_title=title,
        show_heatmap_legend=TRUE,
        column_title_gp=gpar(fontsize=12),
        row_names_gp=gpar(fontsize=2))
dev.off()
