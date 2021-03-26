# The local run part
# visualize the simulation reconstruction results
# the HPC part simu_network_testcluster.R and parallel_job_submit.R need to be finished before this
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(deSolve)
require(igraph)
require(CausalKinetiX)
require(RCy3)
require(ggplot2)
require(reshape2)
# the user will need to modify the path here for local run
comp="/Users/yuewu/"
pardir=paste0(comp,"Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/");
# OR
# pardir=paste0(getwd(),"/")
resdir=paste0(pardir,"result/cluster_simu_test/")
# standard error function
std <- function(x) sd(x)/sqrt(length(x))
# Load the files from HPC
setwd(resdir)
dir_load=paste0(resdir,"parallel_run_big_new/partial/")
storefiles=list.files(path=dir_load,pattern="net\\.CausalKinetiX\\.network\\_\\d+\\.RData")
storefiles %>% str_replace_all(string=.,pattern=fixed("net.CausalKinetiX.network_"),replacement="") %>%
               str_replace_all(string=.,pattern=fixed(".RData"),replacement="") %>%
               as.numeric(.) -> nnetworks_seq
listres=sapply(nnetworks_seq,simplify=FALSE,function(x){
  load(paste0(dir_load,"net.CausalKinetiX.network_",x,".RData"))
  listressep
})
listres=unlist(listres,recursive=FALSE)
# parameters
ncondiseq=c(2,10,40,100)
## network clustering
cytoscapePing()
cytoscapeVersionInfo()
# set the default setting in cytoscape of clustering to "Create new clustered network"=TRUE, "Restore inter-cluster edges after layout"=TRUE
listeval=vector(mode="list")
for(inetwork in nnetworks_seq){#different random network
  for(ncondi in ncondiseq){##different conditions
    filepath=paste0(dir_load,"net.CausalKinetiX_network",inetwork,"_condlen_",ncondi,".txt")
    edgetab=read.table(filepath,sep="\t",header=TRUE)
    edgetab[,"source"]=as.character(edgetab[,"source"])
    edgetab[,"target"]=as.character(edgetab[,"target"])
    createNetworkFromDataFrames(edge=edgetab,title="causalkinetix_network",collection="causalkinetix_network")
    clustermaker_url<-paste0("cluster glay undirectedEdges=true createGroups=false")
    commandsGET(clustermaker_url)
    clusttab=getTableColumns(table="node")
    subtab=clusttab[,c("name","__glayCluster")]
    clusters=unique(subtab[,"__glayCluster"])
    deleteAllNetworks()
    clust_real=listres[[paste0(inetwork,"_",ncondi)]][["clustseq"]]
    subtab[,"name"]=as.numeric(subtab[,"name"])
    clust_est=subtab[order(subtab[,"name"]),"__glayCluster"]
    # recalculate observed real clusters
    clust_real_obs=clust_real[sort(subtab[,"name"])]
    # for each real cluster, find the overlap size of the estimated cluster
    tab_clust=c()
    for(clust in setdiff(unique(clust_real_obs),0)){
      mem_real=which(clust_real_obs==clust)
      overlaped_clusters=table(clust_est[mem_real])
      big_overlap_size=sort(overlaped_clusters,decreasing=TRUE)[1]
      est_clust_size=length(which(clust_est==as.numeric(names(big_overlap_size))))
      exp_match_prop=length(mem_real)/length(clust_real_obs)
      tab_clust=rbind(tab_clust,c(clust,as.numeric(names(big_overlap_size)),big_overlap_size/est_clust_size,exp_match_prop))#length(mem_real)
    }
    colnames(tab_clust)=c("cluster_real","cluster_match","overlap_prop","exp_prop")
    listeval[[paste0(inetwork,"_",ncondi)]]=tab_clust
  }
}
save(listeval,file=paste0("./","net.CausalKinetiX.all.clusterres.bigpartial.new.RData"))
# plot proportion of overlap
clust_overlap_list=list(clust=c(),condi=c(),val=c())
for(ncondi in ncondiseq){##different conditions
  for(inetwork in nnetworks_seq){#different random network
    tab=listeval[[paste0(inetwork,"_",ncondi)]]
    clust_overlap_list$clust=c(clust_overlap_list$clust,tab[,"cluster_real"])
    clust_overlap_list$condi=c(clust_overlap_list$condi,rep(ncondi,times=dim(tab)[1]))
    clust_overlap_list$val=c(clust_overlap_list$val,tab[,"overlap_prop"])
  }
}
clust_overlap_frame=as.data.frame(clust_overlap_list)
clust_overlap_frame$condi=as.factor(clust_overlap_frame$condi)
clust_overlap_frame$clust=as.factor(clust_overlap_frame$clust)
# plot with 2sd
meantab=aggregate(val~clust+condi,data=clust_overlap_frame,mean)
colnames(meantab)=c("clust","condi","mean")
sdtab=aggregate(val~clust+condi,data=clust_overlap_frame,std)
colnames(sdtab)=c("clust","condi","se")
summarytab=merge(meantab,sdtab,by=c("clust","condi"),all=TRUE)
real_ratio=as.data.frame(listeval[[1]][,c("cluster_real","exp_prop")])
colnames(real_ratio)=c("clust","realratio")
summarytab=merge(summarytab,real_ratio,by=c("clust"),all=TRUE)
p<-ggplot(data=summarytab,aes(x=condi,y=mean,color=clust,group=clust))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=mean-se*2,ymax=mean+se*2))+
   geom_hline(aes(yintercept=realratio),linetype="dashed",color="red",size=0.4)+
   facet_wrap(~clust,nrow=1,scales="free_y")+
   labs(title=expression(paste("match ratio",)),x="condition number",y="match ratio")+
   theme(axis.text=element_text(size=10))
ggsave(paste0("partial_overlapratio_boxplot.pdf"),plot=p)
# plot with raw data
clust_overlap_frame=merge(clust_overlap_frame,real_ratio,by=c("clust"),all=TRUE)
p<-ggplot(data=clust_overlap_frame,aes(x=condi,y=val,color=clust,group=clust))+
  geom_point()+
  geom_hline(aes(yintercept=realratio),linetype="dashed",color="red",size=0.4)+
  facet_wrap(~clust,nrow=1,scales="free_y")+
  labs(title=expression(paste("match ratio",)),x="condition number",y="match ratio")+
  theme(axis.text=element_text(size=10))
ggsave(paste0("partial_overlapratio_raw.pdf"),plot=p)
# wilcox t test
stattab=c()
for(clust in unique(clust_overlap_frame[,"clust"])){
  for(cond in unique(clust_overlap_frame[,"condi"])){
    loctab=clust_overlap_frame[clust_overlap_frame[,"clust"]==clust&clust_overlap_frame[,"condi"]==cond,]
    pval=wilcox.test(loctab[,"val"],mu=unique(loctab[,"realratio"]),alternative="greater")$p.value
    stattab=rbind(stattab,c(clust,cond,pval))
  }
}

# accuracy and recovery
liststat=vector(mode="list")
for(inetwork in nnetworks_seq){#different random network
  for(ncondi in ncondiseq){##different conditions
    temp_list=listres[[paste0(inetwork,"_",ncondi)]]
    edgetab_real=temp_list[["simu_network"]]
    edgetab_real=edgetab_real[,c("from","to")]
    edgetab_est=temp_list[["nettab"]]
    edgetab_est=edgetab_est[,c("source","target")]
    colnames(edgetab_est)=c("from","to")
    # there will be edges unobervable so should not be counted as ground truth
    obs_nodes=unique(unlist(edgetab_est))
    edgetab_real=edgetab_real[edgetab_real[,"from"]%in%obs_nodes & edgetab_real[,"to"]%in%obs_nodes,]
    corr_count=sum(apply(edgetab_real,1,function(x){
      loc=which((edgetab_est[,1]==x[1]&edgetab_est[,2]==x[2])|(edgetab_est[,1]==x[2]&edgetab_est[,2]==x[1]))
      length(loc)>0
    }))
    liststat[[paste0(inetwork,"_",ncondi)]]=c(corr_count/dim(edgetab_real)[1],corr_count/dim(edgetab_est)[1])
  }
}
# stat
clust_stat_list=list(condi=c(),recall=c(),precision=c())
for(ncondi in ncondiseq){##different conditions
  for(inetwork in nnetworks_seq){#different random network
    tab=liststat[[paste0(inetwork,"_",ncondi)]]
    clust_stat_list$condi=c(clust_stat_list$condi,ncondi)
    clust_stat_list$recall=c(clust_stat_list$recall,tab[1])
    clust_stat_list$precision=c(clust_stat_list$precision,tab[2])
  }
}
clust_stat_frame=as.data.frame(clust_stat_list)
clust_stat_frame$condi=as.factor(clust_stat_frame$condi)
clust_stat_framelong=melt(clust_stat_frame,id=c("condi"))
# plot with 2sd
meantab=aggregate(value~condi+variable,data=clust_stat_framelong,mean)
colnames(meantab)=c("condi","variable","mean")
sdtab=aggregate(value~condi+variable,data=clust_stat_framelong,std)
colnames(sdtab)=c("condi","variable","se")
summarytab=merge(meantab,sdtab,by=c("condi","variable"),all=TRUE)
p<-ggplot(data=summarytab,aes(x=condi,y=mean,color=variable,group=variable))+
  geom_point()+
  geom_errorbar(aes(ymin=mean-se*2,ymax=mean+se*2))+
  facet_wrap(~variable,nrow=1,scales="free_y")+
  labs(title=expression(paste("precision and recall",)),x="condition number",y="value")+
  theme(axis.text=element_text(size=10))
ggsave(paste0("partial_prec_recall.pdf"),plot=p)
