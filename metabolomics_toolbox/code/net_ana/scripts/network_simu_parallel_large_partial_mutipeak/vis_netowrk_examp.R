# sample visualization of randome sampled network and dynamics
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(deSolve)
require(igraph)
require(ggplot2)
require(reshape2)
require(RCy3)
# the user will need to modify the path here for local run
comp="/Users/yuewu/"
pardir=paste0(comp,"Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/");
# OR
# pardir=paste0(getwd(),"/")
resdir=paste0(pardir,"result/cluster_simu_test/")
dir_load=paste0(resdir,"parallel_run_big_new/multipeak/")
setwd(resdir)
#load data
inetwork=1
ncondi=100
load(paste0(dir_load,"net.CausalKinetiX.network_",inetwork,".RData"))
control_points=c("0","1","2","3")
colors=c("#CCCCCC","#FF0033","#00FF00","#0033FF")
style="BioPAX_0_0"#please use the style you like and have
#
locres=listressep[["1_100"]]
edgetab=locres[["simu_network"]][,c("from","to")]
# recover compound based clusters from peak based cluster table
clust_peak_tab=locres[["clustseq"]]
peak_tab=locres[["peakvec"]]
clust_compd_tab=t(sapply(unique(peak_tab),function(x){
  peakind=which(peak_tab==x)
  clu_label=unique(clust_peak_tab[peakind])
  if(length(clu_label)>1){
    stop("cluster information miss")
  }
  c(x,clu_label)
}))
nettab_ori_edge=data.frame(source=edgetab[,"from"],target=edgetab[,"to"])
nettab_ori_edge[,"source"]=as.character(nettab_ori_edge[,"source"])
nettab_ori_edge[,"target"]=as.character(nettab_ori_edge[,"target"])
nettab_ori_node=as.data.frame(clust_compd_tab)
colnames(nettab_ori_node)=c("id","cluster")
nettab_ori_node[,"id"]=as.character(nettab_ori_node[,"id"])
nettab_ori_node[,"cluster"]=as.character(nettab_ori_node[,"cluster"])
cytoscapePing()
cytoscapeVersionInfo()
deleteAllNetworks()
createNetworkFromDataFrames(edges=nettab_ori_edge,nodes=nettab_ori_node,title="random_network",collection="random_network")#
setVisualStyle(style)
setNodeColorMapping("cluster",control_points,colors,mapping.type="discrete",style.name=style)
setNodeSizeDefault(40)
setNodeShapeDefault('ELLIPSE')
setNodeLabelBypass(nettab_ori_node[,"id"],"")
saveSession("random_network_examp.cys")
deleteAllNetworks()
#sample dynamics
# source("/PATH/TO/network_simu.R")
randseed=1
set.seed(randseed)
n_node=100
nrep=3#number of replicates
krange=c(0,1)#range of kinetic parameters
times=seq(from=0,to=5,by=0.2)#the time sequence
repinfor=list()
repinfor[["peakvec"]]=locres[["peakvec"]]
repinfor[["peakfactor"]]=c(0.3,3)
function_str=react_constr(locres[["simu_network"]],n_node)
cond_vec=rep(seq(ncondi),each=nrep)
unicond=unique(cond_vec)
# initial conditions
iniarra=matrix(NA,nrow=n_node,ncol=length(unicond))
for(cond_i in seq(length(unicond))){
  condvec=runif(n_node,min=krange[1],max=krange[2])
  iniarra[,cond_i]=condvec
}
matdata=ode_simu(odefunc=function_str,cond_vec=cond_vec,times=times,iniarra=iniarra,repinfor=repinfor)
node_seq=c(1,8,13)
replic_seq=1:3
condi_seq=1:3
ntime=length(times)
plotlist=list("time"=c(),"value"=c(),"replic"=c(),"condi"=c(),"nodes"=c())
for(nodeele in node_seq){
  for(condele in condi_seq){
    for(replicele in replic_seq){
      rowind=(condele-1)*length(replic_seq)+replicele
      colind=(nodeele-1)*ntime+seq(ntime)
      plotlist[["time"]]=c(plotlist[["time"]],times)
      plotlist[["value"]]=c(plotlist[["value"]],matdata[rowind,colind])
      plotlist[["replic"]]=c(plotlist[["replic"]],rep(replicele,times=ntime))
      plotlist[["condi"]]=c(plotlist[["condi"]],rep(condele,times=ntime))
      plotlist[["nodes"]]=c(plotlist[["nodes"]],rep(nodeele,times=ntime))
    }
  }
}
plottab=as.data.frame(plotlist)
for(faccol in c("replic","condi","nodes")){
  plottab[,faccol]=as.factor(plottab[,faccol])
}
meantab=aggregate(value~condi+nodes+time,data=plottab,mean)
colnames(meantab)=c("condi","nodes","time","mean")
sdtab=aggregate(value~condi+nodes+time,data=plottab,sd)
colnames(sdtab)=c("condi","nodes","time","sd")
summarytab=merge(meantab,sdtab,by=c("condi","nodes","time"),all=TRUE)
p<-ggplot(data=summarytab,aes(x=time,y=mean,color=condi,group=nodes))+
    geom_point()+
    # geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
    facet_wrap(~nodes,nrow=1,scales="free_y")+
    labs(title=expression(paste("sample dynamic curves",)),x="time",y="value")+
    theme(axis.text=element_text(size=10))
ggsave(paste0("exam_node_dynamics.pdf"),plot=p)
