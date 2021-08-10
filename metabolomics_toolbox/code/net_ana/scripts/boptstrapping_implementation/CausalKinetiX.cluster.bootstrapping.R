## use bootstrapping to estimate stablity in network construction from CausalKinetiX
## to run the code here, the user need proper HPC environments, corresponding submit.sh, and the input data
## For running the bootstrapping locally or on other HPCs, please adapt the script correspondingly
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
dir=getwd()
args=commandArgs(trailingOnly=TRUE)
print(args)
if(length(args)<2){
  stop("At least two argument must be supplied for the range of bootstrapping", call.=FALSE)
}
nbootstrpseq=args[1]:args[2]
featurefilt_pval=0.05#filter the feature score to select connected features
sampind=list("1"=c(1,2,3),"2"=c(4,5,6))
para=list(max.preds=TRUE,expsize=2,#expected number of term
          interactions=TRUE,pen.degree=3,smooth.Y=TRUE,
          integrated.model=FALSE,#derivative based objective function
          screening=30,#screening predictor size
          silent=TRUE)
strlist=c("_peakwise")
for(add_str in strlist){
  if(!dir.exists(paste0("res",add_str))){
    dir.create(paste0("res",add_str))
  }
  matres=readMat(paste0(dir,"/CausalKinetiX_input",add_str,".mat"))
  compds=unlist(matres[["namesall"]])
  time=as.numeric(matres[["timevec"]])
  ntime=length(time)
  matdata=matres[["mat.all"]]
  env=as.numeric(matres[["envvec"]])
  ncompd=length(compds)
  for(ibootstrap in nbootstrpseq){
    set.seed(ibootstrap)
    currdir=paste0("res",add_str,"/",ibootstrap)
    dir.create(currdir)
    # resampling with replacement independently for each compound and separately for two conditions
    matdatanew=matdata
    for(compdi in 1:ncompd){
      colind=((compdi-1)*ntime+1):(compdi*ntime)
      resampind=c(sample(sampind[[1]],3,replace=TRUE),sample(sampind[[2]],3,replace=TRUE))
      matdatanew[unlist(sampind),colind]=matdata[resampind,colind]
    }
    # loop though each compound or peaks
    inforlist<-foreach(compdi=seq(length(compds)))%dopar%{
      compd=compds[compdi]
      print(compd)
      ck.fit<-CausalKinetiX(matdatanew,time-min(time),env,compdi,pars=para)
      scores=ck.fit[["variable.scores"]]
      seleind=which(scores<featurefilt_pval)
      selefeature=ck.fit[["ranking"]][seleind]
      list(detailres=ck.fit,compd=compd,compdi=compdi,selefeature=selefeature)
    }
    save(inforlist,file=paste0("./",currdir,"/rawdata.1",add_str,".RData"))
    ## format the data into network link matrix
    assoc_mat=matrix(0,nrow=ncompd,ncol=ncompd)##row Y col X
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
    for(i in 1:ncompd){
      edgeind=which(assoc_mat_sym[i,]>0)
      sourcenod=c(sourcenod,rep(compds[i],times=length(edgeind)))
      targetnode=c(targetnode,compds[edgeind])
      assovec=c(assovec,assoc_mat_sym[i,edgeind])
    }
    nettab=data.frame(source=sourcenod,target=targetnode,association=assovec)
    save(nettab,file=paste0("./",currdir,"/net.CausalKinetiX",add_str,".RData"))
    write.table(nettab,file=paste0("./",currdir,"/net.CausalKinetiX",add_str,".txt"),sep="\t",row.names=FALSE,quote=FALSE)
  }
}
