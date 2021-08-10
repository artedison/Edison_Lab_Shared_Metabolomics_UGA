# test cases listed here

test_that("random network simulation proportion correct", {
  n_node=200
  cluster_size_vec=c(50,40,40)
  p_edge_vec=c(0.6,0.3)
  p_regu=0.2
  p_comb=0.2
  krange=c(0,10)
  randseed=1
  inforlist=network_constr(n_node=n_node,cluster_size_vec=cluster_size_vec,p_edge_vec=p_edge_vec,p_regu=p_regu,p_comb=p_comb,krange=krange,randseed=randseed)
  edgetab=inforlist[["edgetab"]]
  clustseq=inforlist[["clustseq"]]
  #
  exp_size=c(n_node-sum(cluster_size_vec),cluster_size_vec)
  clus_flag=all(sort(exp_size)-sort(table(clustseq))==0)
  if(!clus_flag){
    print("cluster size different")
  }
  p_regu_meas=(table(edgetab[,"type"])/nrow(edgetab))[2]
  regu_p_flag=(p_regu-p_regu_meas)<0.01
  if(!regu_p_flag){
    print("regu p different")
  }
  p_comb_meas=length(which(!is.na(edgetab[,"combination"])))/n_node/2
  comb_p_flag=(p_comb-p_comb_meas)<0.01
  if(!comb_p_flag){
    print("comb p different")
  }
  p_rev_meas=length(which(!is.na(edgetab[,"rev"])))/nrow(edgetab)
  rev_p_flag=(0.33-p_rev_meas)<0.01
  if(!rev_p_flag){
    print("rev p different")
  }
  regudir=table(edgetab[,"regudir"])
  p_reg_meas=regudir[1]/sum(regudir)
  reg_p_flag=(0.5-p_reg_meas)<0.01
  if(!reg_p_flag){
    print("regu direction p different")
  }
  stattest=ks.test(edgetab[,"k"],"punif",0,10)
  unif_flag=stattest$p.value>0.1 && min(edgetab[,"k"])>0 && max(edgetab[,"k"])<10
  if(!unif_flag){
    print("k distribution not uniform")
  }
  expect_true(clus_flag&regu_p_flag&comb_p_flag&rev_p_flag&reg_p_flag&unif_flag)
})

test_that("random network simulation regulation target connection correct", {
  n_node=200
  cluster_size_vec=c(50,40,40)
  p_edge_vec=c(0.6,0.3)
  p_regu=0.2
  p_comb=0.2
  krange=c(0,10)
  randseed=1
  inforlist=network_constr(n_node=n_node,cluster_size_vec=cluster_size_vec,p_edge_vec=p_edge_vec,p_regu=p_regu,p_comb=p_comb,krange=krange,randseed=randseed)
  edgetab=inforlist[["edgetab"]]
  clustseq=inforlist[["clustseq"]]
  regu_ind=which(edgetab[,"type"]=="regu")
  expect_true(all(edgetab[regu_ind,"to"]==edgetab[edgetab[regu_ind,"regutarget"],"to"]))
})

test_that("regulation not reversible", {
  n_node=200
  cluster_size_vec=c(50,40,40)
  p_edge_vec=c(0.6,0.3)
  p_regu=0.2
  p_comb=0.2
  krange=c(0,10)
  randseed=1
  inforlist=network_constr(n_node=n_node,cluster_size_vec=cluster_size_vec,p_edge_vec=p_edge_vec,p_regu=p_regu,p_comb=p_comb,krange=krange,randseed=randseed)
  edgetab=inforlist[["edgetab"]]
  clustseq=inforlist[["clustseq"]]
  regu_ind=which(edgetab[,"type"]=="regu")
  expect_true(all(is.na(edgetab[regu_ind,"rev"])))
})

test_that("random network simulation regulation target connection correct for small network", {
  n_node=30
  cluster_size_vec=c(14,10)
  p_edge_vec=c(0.2,0.05)
  p_regu=0.2
  p_comb=0.2
  krange=c(0,1)
  randseed=1
  inforlist=network_constr(n_node=n_node,cluster_size_vec=cluster_size_vec,p_edge_vec=p_edge_vec,p_regu=p_regu,p_comb=p_comb,krange=krange,randseed=randseed)
  edgetab=inforlist[["edgetab"]]
  clustseq=inforlist[["clustseq"]]
  regu_ind=which(edgetab[,"type"]=="regu")
  expect_true(all(edgetab[regu_ind,"to"]==edgetab[edgetab[regu_ind,"regutarget"],"to"]))
})

test_that("random network simulation combined reaction sharing qualities", {
  n_node=200
  cluster_size_vec=c(50,40,40)
  p_edge_vec=c(0.6,0.3)
  p_regu=0.2
  p_comb=0.2
  krange=c(0,10)
  randseed=1
  inforlist=network_constr(n_node=n_node,cluster_size_vec=cluster_size_vec,p_edge_vec=p_edge_vec,p_regu=p_regu,p_comb=p_comb,krange=krange,randseed=randseed)
  edgetab=inforlist[["edgetab"]]
  clustseq=inforlist[["clustseq"]]
  comb_inds=which(!is.na(edgetab[,"combination"]))
  flag_all=TRUE
  for(ind in comb_inds){
    ind2=edgetab[ind,"combination"]
    kflag=edgetab[ind2,"k"]==edgetab[ind,"k"]
    toflag=edgetab[ind2,"to"]==edgetab[ind,"to"]
    typeflag=edgetab[ind2,"type"]=="reac"&&edgetab[ind,"type"]=="reac"
    combflag=edgetab[ind2,"combination"]==ind
    revflag=is.na(edgetab[ind2,"rev"])&&is.na(edgetab[ind,"rev"])
    if(!(kflag&&toflag&&typeflag&&combflag&&revflag)){
      flag_all=FALSE
      break
    }
  }
  expect_true(flag_all)
})

test_that("generating equation from network small test", {
  n_node=8
  from=c(1,3,2,3,4,5,6,7)
  to=c(2,2,4,4,5,1,4,4)
  type=rep("reac",times=8)
  type[4]="regu"
  combination=rep(NA,times=8)
  combination[c(3,7)]=c(7,3)
  rev=rep(NA,times=8)
  rev[2]=1
  regudir=rep(NA,times=8)
  regudir[4]=1
  regutarget=rep(NA,times=8)
  regutarget[4]=3
  k=c(1,2,3,2,1,2,3,1)
  edgetab=data.frame(from=from,to=to,type=type,combination=combination,rev=rev,regudir=regudir,regutarget=regutarget,k=k)
  function_str=react_constr(edgetab,n_node)
  function_str_exp="odef<-function(t,y,parms) {\ndy=numeric(8)\ndy[1]=-1*y[1]+1*2*y[5]\ndy[2]=+1*1*y[1]+1*2*y[3]-1*2*y[2]-3*y[3]^1*y[2]*y[6]\ndy[3]=-2*y[3]+2*y[2]\ndy[4]=+2*3*y[3]^1*y[2]*y[6]-1*y[4]+1*1*y[7]\ndy[5]=+1*1*y[4]-2*y[5]\ndy[6]=-3*y[3]^1*y[2]*y[6]\ndy[7]=-1*y[7]\ndy[8]=0\nreturn(list(dy))\n}"
  expect_identical(function_str,function_str_exp)
})

test_that("ode mass conserve", {
  set.seed(1)
  n_node=8
  from=c(1,3,2,3,4,5,6,7)
  to=c(2,2,4,4,5,1,4,4)
  type=rep("reac",times=8)
  type[4]="regu"
  combination=rep(NA,times=8)
  combination[c(3,7)]=c(7,3)
  rev=rep(NA,times=8)
  rev[2]=1
  regudir=rep(NA,times=8)
  regudir[4]=1
  regutarget=rep(NA,times=8)
  regutarget[4]=3
  k=c(1,2,3,2,1,2,3,1)
  edgetab=data.frame(from=from,to=to,type=type,combination=combination,rev=rev,regudir=regudir,regutarget=regutarget,k=k)
  odefunc=react_constr(edgetab,n_node)
  #
  cond_vec=c(1,1,1,2,2,2,3,3,3)
  krange=c(0,10)
  times=seq(from=0,to=1,by=0.01)
  unicond=unique(cond_vec)
  iniarra=matrix(NA,nrow=n_node,ncol=length(unicond))
  for(cond_i in seq(length(unicond))){
    condvec=runif(n_node,min=krange[1],max=krange[2])
    iniarra[,cond_i]=condvec
  }
  matdata=ode_simu(odefunc=odefunc,cond_vec=cond_vec,times=times,iniarra=iniarra)
  #test the conserve for each condition&replicate&time
  mass_cons_thre=0.01
  array_size=dim(matdata)
  ntime=length(times)
  mass=NA;
  flag_all=TRUE;
  for(repi in 1:array_size[1]){
    for(timei in 1:ntime){
      timeslic=seq(from=timei,to=array_size[2],by=ntime)
      mass_t=sum(matdata[repi,timeslic])
      if(timei==1){
        mass=mass_t
      }else{
        if(abs(mass-mass_t)/mass>mass_cons_thre){
          warning(paste0("mass not conserved at replicate ",repi," time ",timei));
          flag_all=FALSE;
          break;
        }
      }
    }
  }
  expect_true(flag_all)
})

test_that("ode duplicated peak simu", {
  set.seed(1)
  n_node=8
  from=c(1,3,2,3,4,5,6,7)
  to=c(2,2,4,4,5,1,4,4)
  type=rep("reac",times=8)
  type[4]="regu"
  combination=rep(NA,times=8)
  combination[c(3,7)]=c(7,3)
  rev=rep(NA,times=8)
  rev[2]=1
  regudir=rep(NA,times=8)
  regudir[4]=1
  regutarget=rep(NA,times=8)
  regutarget[4]=3
  k=c(1,2,3,2,1,2,3,1)
  edgetab=data.frame(from=from,to=to,type=type,combination=combination,rev=rev,regudir=regudir,regutarget=regutarget,k=k)
  odefunc=react_constr(edgetab,n_node)
  #
  cond_vec=c(1,1,1,2,2,2,3,3,3)
  krange=c(0,10)
  times=seq(from=0,to=1,by=0.01)
  unicond=unique(cond_vec)
  iniarra=matrix(NA,nrow=n_node,ncol=length(unicond))
  for(cond_i in seq(length(unicond))){
    condvec=runif(n_node,min=krange[1],max=krange[2])
    iniarra[,cond_i]=condvec
  }
  # peak repeat information
  nrepeat=sample(c(1,2,3,4,5),n_node,replace=TRUE)
  repinfor=list()
  repinfor[["peakvec"]]=rep(seq(n_node),times=nrepeat)
  repinfor[["peakfactor"]]=c(0.3,3)
  matdata=ode_simu(odefunc=odefunc,cond_vec=cond_vec,times=times,iniarra=iniarra,repinfor=repinfor)
  # dimension of matdata
  flag_dim=(dim(matdata)[2]/length(times)-length(repinfor[["peakvec"]]))==0
  # correlation of peaks in the same group
  matdata_refor=c(matdata)
  dim(matdata_refor)=c(9*101,20)
  cormat=cor(matdata_refor)
  cormat_thre=cormat>0.9
  counts=colSums(cormat_thre)
  flag_group=all(sapply(unique(repinfor[["peakvec"]]),function(x){
    ind=which(x==repinfor[["peakvec"]])
    all(length(ind)==counts[ind])
  }))
  expect_true(flag_dim&flag_group)
})
