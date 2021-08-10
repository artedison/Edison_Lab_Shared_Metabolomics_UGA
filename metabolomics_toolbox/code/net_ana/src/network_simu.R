#' Random simulate a network with clusters.
#'
# The network will contain dense clusters. Nodes represent compounds and edges represent reactions or regulations. Both nodes and edges will be randomly generated. Direction and types of reactions will also be random. For each nodes, one pair of random selected edges of same direction will be selected for combination. Regulation target will be random selected (50% positive|negative regulations). Kinetic parameter will be random generated.
#'
#' @param n_node int. number of nodes to be included in the network. must be provided
#' @param cluster_size_vec int array. length(cluster_size_vec)=nclusters. size of each clusters. must be provided
#' @param p_edge_vec float array length(p_edge_vec)==2. the probablity of edge connection in and out of clusters. c(p_in,p_out). must be provided
#' @param p_regu float. probablity of a edge to be regulation edge. default 0.2
#' @param p_comb float. probablity to include combinations (2 reactants) for a product. default 0.2
#' @param krange float array. length(krange)==2. the parameter range for k. default c(0,10)
#' @param randseed int. the ranomd number seed. default 1.
#' @return a list containing:
#' edgetab. The network edge table. Columns: from node, to node, node type ('reac','regu'),  combination (the reaction ind that combined with this reaction, will share k), rev (reversibility, 1 if reversible), regulation directions (1 up-regulation, -1 down-regulation), regulation target (the reaction that regulated, will share the same to), kinetic parameter. Correspinding column names: from, to, type, combination, rev, regudir, regutarget, k. nrow(edgetab)=nedges
#' clustseq. the vector indicating which cluster the nodes belongs
#' @details
#'  # Caveats
#'  k for regulation edges doesn't contribute to the equation
#'  reversible reactions share the same k
#'  combined reaction are irreversible
#'  No combination of regulation and reversible reactions
#'  Author: YUE WU (yue.wu@uga.edu)
#'
network_constr<-function(n_node=NULL,cluster_size_vec=NULL,p_edge_vec=NULL,p_regu=0.2,p_comb=0.2,krange=c(0,10),randseed=1){
  if(is.null(n_node)){
    stop("please provide the number of nodes")
  }
  if(is.null(cluster_size_vec)){
    stop("please provide the cluster size vector")
  }
  if(is.null(p_edge_vec)){
    stop("please provide the probablity of edge connection vector (in and outside of clusters)")
  }
  set.seed(randseed)
  nodeseq=seq(from=1,to=n_node,by=1)
  clustseq=rep(0,times=n_node)
  # sampling nodes to clusters
  cluster_labels=seq(length(cluster_size_vec))
  for(clus_size_i in cluster_labels){
    clus_size=cluster_size_vec[clus_size_i]
    ind_incluster=sample(nodeseq[clustseq==0],clus_size,replace=FALSE)
    clustseq[ind_incluster]=clus_size_i
  }
  # sampling edges within clusters (including all nodes)
  # all nodes are first connected with a small p and then nodes in each cluster are connected with a higher p
  cluster_size_vec_all=c(n_node,cluster_size_vec)
  cluster_label_all=c(0,cluster_labels)
  edgetab=c()
  for(clus_size_i in seq(length(cluster_label_all))){
    cluster_size=cluster_size_vec_all[clus_size_i]
    cluster_label=cluster_label_all[clus_size_i]
    if(cluster_label==0){
      p_edge=p_edge_vec[2]
      alledges=combn(nodeseq,2)
    }else{
      p_edge=p_edge_vec[1]
      alledges=combn(nodeseq[clustseq==cluster_label],2)
    }
    nalledges=ncol(alledges)#max number of edge
    nedges=floor(p_edge*nalledges)
    temptab=alledges[,sample(seq(nalledges),nedges)]
    edgetab=rbind(edgetab,t(temptab))
  }
  # remove duplicated edges
  remmask=duplicated(as.data.frame(t(apply(edgetab,1,function(x){
    sort(x)
  }))))
  edgetab=edgetab[!remmask,]
  edgetab=as.data.frame(edgetab)
  colnames(edgetab)=c("from","to")
  # sampling edge types
  nedge_all=nrow(edgetab)
  edgetab$type=sample(c("reac","regu"),nedge_all,replace=TRUE,prob=c(1-p_regu,p_regu))
  # sampling edgedirections
  flip_edge=which(sample(c(1,0),nedge_all,replace=TRUE)==1)
  temp=edgetab[flip_edge,"from"]
  edgetab[flip_edge,"from"]=edgetab[flip_edge,"to"]
  edgetab[flip_edge,"to"]=temp
  # sampling combinations
  node_edge_seq_ind=rep(0,times=n_node)
  node_edge_seq_ind[unique(c(edgetab[,"from"],edgetab[,"to"]))]=1#nodes with at least one edge
  node_edge_seq=nodeseq[which(node_edge_seq_ind==1)]
  sele_nodes=sample(node_edge_seq,floor(p_comb*length(node_edge_seq)),replace=FALSE)
  combvec=rep(NA,times=nedge_all)
  kvec_comb=rep(NA,times=nedge_all)
  for(node in sele_nodes){
    in_reac_edge=which(edgetab[,"to"]==node&edgetab[,"type"]=="reac")
    if(length(in_reac_edge)<2){
      next
    }
    pairs=combn(in_reac_edge,2)
    combpair=pairs[,sample(seq(ncol(pairs)),1)]
    combvec[combpair[1]]=combpair[2]
    combvec[combpair[2]]=combpair[1]
    randk=runif(1,min=krange[1],max=krange[2])
    kvec_comb[c(combpair[1],combpair[2])]=randk
  }
  edgetab$combination=combvec
  # add reversible reactions
  p_rev=1/3
  rev_vec=rep(NA,times=nedge_all)
  rev_vec[sample(seq(nedge_all)[is.na(combvec)&edgetab$type=="reac"],floor(nedge_all*p_rev),replace=FALSE)]=1#combined reaction are irresible
  edgetab$rev=rev_vec
  # random generate direction&target of regulations
  regu_direc_vec=rep(NA,times=nedge_all)
  regu_target_vec=rep(NA,times=nedge_all)
  indregu=which(edgetab$type=="regu")
  indreac=which(edgetab$type=="reac")
  regu_direc_vec[indregu]=sample(c(-1,1),length(indregu),replace=TRUE)
  for(singregu_ind in indregu){
    targetnode=edgetab[singregu_ind,"to"]
    tar_reactions=which(edgetab[,"to"]==targetnode&edgetab[,"type"]=="reac"&is.na(edgetab[,"rev"]))#no regulation + reversible
    if(length(tar_reactions)>0){
      if(length(tar_reactions)==1){
        regu_target_vec[singregu_ind]=tar_reactions
      }else{
        regu_target_vec[singregu_ind]=sample(tar_reactions,1,replace=TRUE)
      }
    }
  }
  naflag=is.na(regu_direc_vec)|is.na(regu_target_vec)
  regu_direc_vec[naflag]=NA
  regu_target_vec[naflag]=NA
  edgetab$regudir=regu_direc_vec
  edgetab$regutarget=regu_target_vec
  edgetab[naflag,"type"]="reac"
  # random generate reaction parameter k
  kvec=runif(nedge_all,min=krange[1],max=krange[2])
  indcomb=which(!is.na(kvec_comb))
  kvec[indcomb]=kvec_comb[indcomb]
  edgetab$k=kvec
  return(list(edgetab=edgetab,clustseq=clustseq))
}

#' Constructing ODE set based on the network edge table.
#'
#' This function will take the network edge table and fomulate the ODE equations as a function. Reactions will be generated for the two nodes of each edges. The ODE function will be automatically generated and returned. This function should be determinstic.
#'
#' @param edgetab dataframe. the network edge table from [network_constr()]. Column names: from, to, type, combination, rev, regudir, regutarget, k. nrow(edgetab)=nedges. must be provided
#' @param n_node int. number of nodes to be included in the network. must be provided
#'
#' @return function_str. the ODE function string to be included for ODE simulation
#' @details
#'  Author: YUE WU (yue.wu@uga.edu)
#'
react_constr<-function(edgetab=NULL,n_node=NULL){
  if(is.null(edgetab)){
    stop("please provide the edge table")
  }
  if(is.null(n_node)){
    stop("please provide the number of nodes")
  }
  args=paste0(c("t","y","parms"),collapse=",")
  # initialize the function equation array
  reactants=paste0("y[",seq(n_node),"]",sep="")
  reactions=paste0("d",reactants,"=",sep="")
  # update k from regulations
  k_upd=paste0(edgetab[,"k"])
  for(reguind in which(edgetab[,"type"]=="regu")){
    locinfor=edgetab[reguind,]
    target=locinfor[,"regutarget"]
    k_upd[target]=paste0(k_upd[target],"*",reactants[locinfor[,"from"]],"^",locinfor[,"regudir"])
    if(!is.na(edgetab[target,"combination"])){
      target_comb=locinfor[target,"combination"]
      k_upd[target_comb]=paste0(k_upd[target_comb],"*",reactants[locinfor[,"from"]],"^",locinfor[,"regudir"])
    }
  }
  # contruct reaction from each row
  # keep only one combination part
  comb_vec=rep(NA,times=nrow(edgetab))
  for(reacind in which(edgetab[,"type"]=="reac")){
    if(!is.na(comb_vec[reacind])){
      next
    }
    locinfor=edgetab[reacind,]
    from=locinfor[,"from"]
    to=locinfor[,"to"]
    tofactor=1
    if(!is.na(locinfor[,"combination"])){
      pairedge=locinfor[,"combination"]
      from=c(from,edgetab[pairedge,"from"])
      comb_vec[pairedge]=1
      tofactor=2
    }
    reaction_term=paste0(k_upd[reacind],"*",paste0(reactants[from],collapse="*"))
    reactions[to]=paste0(reactions[to],"+",tofactor,"*",reaction_term)
    reactions[from]=paste0(reactions[from],"-",reaction_term)
    if(!is.na(locinfor[,"rev"])){
      reaction_term=paste0(k_upd[reacind],"*",reactants[to])
      reactions[from]=paste0(reactions[from],"+",reaction_term)
      reactions[to]=paste0(reactions[to],"-",tofactor,"*",reaction_term)
    }
  }
  # considering singular nodes and give 0 derivatives
  singu_nodes_ind=str_which(string=reactions,pattern="=$")
  reactions[singu_nodes_ind]=paste0(reactions[singu_nodes_ind],"0")
  body=paste0(c(paste0("dy=numeric(",n_node,")"),
                reactions,
                "return(list(dy))"),collapse="\n")
  function_str=paste('odef<-function(',args,') {\n',body,'\n}',sep='')
  return(function_str)
}
#' Simulate the ODE set of the network.
#'
#' The ODE function will be used to simulate dynamics in the network under different conditions and add reasonable noises
#'
#' @param odefunc string. the ode function to be used. must be provided
#' @param cond_vec int array. length(cond_vec)=ncondition*nreplicate. the condition vector. Simulating different initial conditions representing different experimental conditions. e.g. c(1,1,1,2,2,2,3,3,3) three conditions and each with three replicates. must be provided.
#' @param times float array. length(times)=ntimes. the time vector. must be provided.
#' @param iniarra float matrix. dim(iniarra)=c(nnodes,nconditions). the initial condtion matrix for each condtion (column) and each nodes (row). must be provided
#' @param repinfor list. the information used for simulating repeat peaks of the same species. default NULL.
#'            peakvec: int array. the repeat index of each species.
#'            peakfactor: int array. length(peakfactor)=2. The range of factors for value peak/species.
#'
#' @return matdata the ODE simulation result stored as a matrix (ncondition*(nnodes*ntime))
#' @details
#'  Author: YUE WU (yue.wu@uga.edu)
#'
ode_simu<-function(odefunc,cond_vec,times,iniarra,repinfor=NULL){
  if(is.null(odefunc)){
    stop("please provide the ode function")
  }
  if(is.null(cond_vec)){
    stop("please provide the condition vector")
  }
  if(is.null(times)){
    stop("please provide the time vector")
  }
  if(is.null(iniarra)){
    stop("please provide the intial condition vector")
  }
  eval(parse(text=odefunc))
  uniqcond=unique(cond_vec)
  matdata=c()
  if(!is.null(repinfor)){
    peakind=repinfor[["peakvec"]]
    facrange=repinfor[["peakfactor"]]
    factorvec=unlist(sapply(unique(peakind),function(x){
      ind_peak_group=which(peakind==x)
      rat_fact=runif(length(ind_peak_group),min=facrange[1],max=facrange[2])
    }))
  }
  for(condi in seq(length(uniqcond))){
    cond=uniqcond[condi]
    odesolu=ode(y=iniarra[,condi],times=times,func=odef,parms=NULL)
    odesolu=odesolu[,2:dim(odesolu)[2]]
    if(!is.null(repinfor)){
      odesolu_list=sapply(unique(peakind),function(x){
        ind_peak_group=which(peakind==x)
        tempmat=sapply(factorvec[ind_peak_group],function(y){
          odesolu[,x]*y
        })
      })
      odesolu=Reduce('cbind',odesolu_list)
    }
    sd_ep=0.02*apply(odesolu,2,sd)+10^(-7)
    # sd_ep=rep(0,times=dim(odesolu)[2])
    for(repi in which(cond_vec==cond)){
      sdmat=sapply(sd_ep,function(x){
        rnorm(length(times),mean=0,sd=x)
      })
      odesolu_noise=odesolu+sdmat
      solutionvec=c(odesolu_noise)
      matdata=rbind(matdata,solutionvec)
    }
  }
  return(matdata)
}
