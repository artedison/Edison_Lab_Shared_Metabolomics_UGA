# script for use simulated network to test/evaluate CausalKinetiX and clustering
# the simulation will produce dynamics based random networks with clusters
# CausalKinetiX is then used to recover the edges and clusters in the network
## to run the code here, the user need proper HPC environments, corresponding submit.sh, and the input data
## For running the bootstrapping locally or on other HPCs, please adapt the script correspondingly
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(deSolve)
require(igraph)
require(CausalKinetiX)
require(ggplot2)
require(reshape2)
require(foreach)
require(doMC)
registerDoMC(cores=20)
source("/PATH/TO/network_simu.R")#the corrresponding functions on network reconstruction
#
dir=getwd()
args=commandArgs(trailingOnly=TRUE)
print(args)
if(length(args)<1){
  stop("At least one argument must be supplied for the network number", call.=FALSE)
}
#
inetwork=as.numeric(args)#different random network controled within different parallel scripts
n_node=100#number of nodes
compds=seq(n_node)
cluster_size_vec=c(40,20,20)#cluter size
p_edge_vec=c(0.15,0.015)#probability for within cluter edges and outside cluster edges
p_regu=0.2#probability for regulation
p_comb=0.2#probability for reaction combination
krange=c(0,1)#range of kinetic parameters
ncondiseq=c(100,40,10,2)#number of conditions in simulation
nrep=3#number of replicates
times=seq(from=0,to=5,by=0.2)#the time sequence
ntimes=length(times)
randseed=1
set.seed(randseed)
# partial observation of nodes. the user can change the proportion of observation here or even make it completely observable
nobs=floor(n_node/2)#number of observed nodes
obs_seq=sample(n_node,nobs,replace=FALSE)
#
para=list(max.preds=TRUE,expsize=2,#expected number of term in the summation
          interactions=TRUE,pen.degree=3,smooth.Y=TRUE,
          integrated.model=FALSE,#derivative based objective function
          screening=30,#screening predictor size
          silent=TRUE)
featurefilt_pval=0.05#filter the feature score to select connected features
#
listressep=vector(mode="list")
for(ncondi in ncondiseq){##different conditions
  templist=vector(mode="list")
  inforlist=network_constr(n_node=n_node,cluster_size_vec=cluster_size_vec,p_edge_vec=p_edge_vec,p_regu=p_regu,p_comb=p_comb,krange=krange,randseed=randseed+inetwork)
  edgetab=inforlist[["edgetab"]]
  templist[["clustseq"]]=inforlist[["clustseq"]]
  templist[["simu_network"]]=edgetab
  # df.g=graph.data.frame(d=edgetab[,c("from","to")],directed=TRUE)
  # plot(df.g,vertex.label=V(df.g)$name)
  function_str=react_constr(edgetab,n_node)
  cond_vec=rep(seq(ncondi),each=nrep)
  unicond=unique(cond_vec)
  # initial conditions
  iniarra=matrix(NA,nrow=n_node,ncol=length(unicond))
  for(cond_i in seq(length(unicond))){
    condvec=runif(n_node,min=krange[1],max=krange[2])
    iniarra[,cond_i]=condvec
  }
  matdata=ode_simu(odefunc=function_str,cond_vec=cond_vec,times=times,iniarra=iniarra)
  obs_ind=unlist(sapply(obs_seq,simplify=FALSE,function(x){
    seq(from=ntimes*(x-1)+1,to=ntimes*x)
  }))
  matdata_obs=matdata[,obs_ind]
  
  # searching for edges
  inforlist<-foreach(compdi=seq(nobs))%dopar%{#
    compd=obs_seq[compdi]
    cat(compd)
    ck.fit<-CausalKinetiX(matdata_obs,times-min(times),cond_vec,compdi,pars=para)#
    scores=ck.fit[["variable.scores"]]
    seleind=which(scores<featurefilt_pval)
    selefeature=ck.fit[["ranking"]][seleind]
    list(detailres=ck.fit,compd=compd,compdi=compdi,selefeature=selefeature)
  }
  save(inforlist,file=paste0("./","storedata_network_",inetwork,"_condlen_",ncondi,".RData"))
  # templist[["causalkinetix_res"]]=inforlist
  
  ## format the data into network link matrix
  assoc_mat=matrix(0,nrow=nobs,ncol=nobs)##row Y col X
  resstorelist=vector(mode="list")
  for(resi in seq(length(inforlist))){
    tempifor=inforlist[[resi]]
    assoc_mat[tempifor[["compdi"]],tempifor[["selefeature"]]]=1
    resstorelist[[tempifor[["compdi"]]]]=tempifor[["detailres"]]
  }
  # scale the histogram
  assoc_mat_sym=(assoc_mat+t(assoc_mat))/2
  assoc_mat_sym[lower.tri(assoc_mat_sym)]=0;
  sourcenod=c()
  targetnode=c()
  assovec=c()
  for(i in 1:nobs){
    edgeind=which(assoc_mat_sym[i,]>0)
    sourcenod=c(sourcenod,rep(obs_seq[i],times=length(edgeind)))
    targetnode=c(targetnode,obs_seq[edgeind])
    assovec=c(assovec,assoc_mat_sym[i,edgeind])
  }
  nettab=data.frame(source=sourcenod,target=targetnode,association=assovec)
  save(nettab,file=paste0("./","net.CausalKinetiX_network_",inetwork,"_condlen_",ncondi,".RData"))
  write.table(nettab,file=paste0("./","net.CausalKinetiX_network",inetwork,"_condlen_",ncondi,".txt"),sep="\t",row.names=FALSE,quote=FALSE)
  templist[["nettab"]]=nettab
  listressep[[paste0(inetwork,"_",ncondi)]]=templist
}
save(listressep,file=paste0("./","net.CausalKinetiX.network_",inetwork,".RData"))
