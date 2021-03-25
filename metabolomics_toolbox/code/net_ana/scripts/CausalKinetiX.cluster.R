# this script will construct the main network by CausalKinetiX based on time series Metabolic measurement
# running the CausalKinetiX block will take multiple core and will take long time.
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(R.matlab)
require(ggplot2)
require(dplyr)
require(CausalKinetiX)
require(foreach)
require(doMC)
registerDoMC(cores=10)
# the user will need to modify the path here for local run
comp="/Users/yuewu/"
pardir=paste0(comp,"Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/");
# OR
# pardir=paste0(getwd(),"/")
datadir=paste0(pardir,"result_data/")
resdir=paste0(pardir,"result/function_cluster/")
#
set.seed(1)
setwd(resdir)
add_str="peakwise"#""
matres=readMat(paste0(datadir,"CausalKinetiX_input_peakwise.mat"))
compds=unlist(matres[["namesall"]])
matdata=matres[["mat.all"]]
time=as.numeric(matres[["timevec"]])
env=as.numeric(matres[["envvec"]])
para=list(max.preds=TRUE,expsize=2,#expected number of term in the summation
          interactions=TRUE,pen.degree=3,smooth.Y=TRUE,
          integrated.model=FALSE,#derivative based objective function
          screening=30,#screening predictor size
          silent=TRUE)
ncompd=length(compds)
featurefilt_pval=0.05#filter the feature score to select connected features

# loop though each compound or peaks (as Y)
inforlist<-foreach(compdi=seq(length(compds)))%dopar%{#
  compd=compds[compdi]
  cat(compd)
  ck.fit<-CausalKinetiX(matdata,time-min(time),env,compdi,pars=para)#
  scores=ck.fit[["variable.scores"]]
  seleind=which(scores<featurefilt_pval)
  selefeature=ck.fit[["ranking"]][seleind]
  list(detailres=ck.fit,compd=compd,compdi=compdi,selefeature=selefeature)
}
save(inforlist,file=paste0("rawdata.1",add_str,".RData"))
## format the data into network link matrix
assoc_mat=matrix(0,nrow=ncompd,ncol=ncompd)##row Y col X
resstorelist=vector(mode="list")
for(resi in seq(length(inforlist))){
  tempifor=inforlist[[resi]]
  assoc_mat[tempifor[["compdi"]],tempifor[["selefeature"]]]=1
  resstorelist[[tempifor[["compdi"]]]]=tempifor[["detailres"]]
}
# symmetric the matrix and keep the up-triagle
assoc_mat_sym=(assoc_mat+t(assoc_mat))/2
assoc_mat_sym[lower.tri(assoc_mat_sym)]=0;
sourcenod=c()
targetnode=c()
assovec=c()
for(i in 1:ncompd){
  edgeind=which(assoc_mat_sym[i,]>0)
  sourcenod=c(sourcenod,rep(compds[i],times=length(edgeind)))
  targetnode=c(targetnode,compds[edgeind])
  assovec=c(assovec,assoc_mat_sym[i,edgeind])
}
nettab=data.frame(source=sourcenod,target=targetnode,association=assovec)
save(nettab,file=paste0("net.CausalKinetiX",add_str,".RData"))
write.table(nettab,file=paste0("net.CausalKinetiX",add_str,".txt"),sep="\t",row.names=FALSE,quote=FALSE)
