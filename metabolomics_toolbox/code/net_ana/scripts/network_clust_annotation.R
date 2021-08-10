# visualize spectra cluster and map to GISSMO database
# The codes here depends on R using Cytoscape through RCy3, R using Python through reticulate, R using MATLAB through R.matlab, and API to GISSMO.
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(R.matlab)
require(RCy3)
require(httr)
require(jsonlite)
require(reticulate)

# the user will need to modify the path here for local run
comp="/Users/yuewu/"
pardir=paste0(comp,"Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/");
# OR
# pardir=paste0(getwd(),"/")
datadir=paste0(pardir,"result_data/")
resdir=paste0(pardir,"result/annotation/mapping/")
setwd(resdir)
Sys.setenv(PATH=paste(Sys.getenv("PATH"),"/Applications/MATLAB_R2018b.app/bin/",sep=":"))

# control parameters
mini_clust_size=3#minimum size of cluster to visualize
mini_match_prop=0.5#minimum proportion of mapped peaks to visualize. one of the OR condition
mini_match_num=5#minimum number of mapped peaks to visualize. one of the OR condition
show_samp_i=1#the visualization will be on sample one.
maxmatch=5#visualize the first 5 peaks

# control cytoscape to obtain clusters
cytoscapePing()
cytoscapeVersionInfo()
file.copy(paste0(datadir,"corr_networkcluster.cys"),paste0(resdir,"corr_networkcluster.copy.cys"))
openSession(paste0(resdir,"corr_networkcluster.copy.cys"))

setCurrentView("causalkinetix_corr.txt--clustered_mcl5")
Sys.sleep(2)#an abstract time of wait to make sure the changes has been made in Cytoscape
nodetab=getTableColumns(table="node")
closeSession(FALSE)

# search peaks on gissmo through API
# conda_create("clust_network_match")
use_condaenv(condaenv="clust_network_match",required=TRUE)
# conda_install("clust_network_match","requests")
### API parameters
peaksearch_url="http://gissmo.nmrfam.wisc.edu/peak_search"
apipara=dict('peak_type'='standard',# Options are 'standard' or 'GSD' for deconvoluted/GSD picking
             'threshold'=0.01,# matching threshold (ppm)
             'frequency'=as.integer(600),# MHz
             'rs'="",#Specify resonances as space separated text string. will be updated later
             'json'=TRUE)

# filter of clusters
clust_vec=nodetab[,"__mclCluster"]
clust_count=table(clust_vec)
clust_filter=names(clust_count)[clust_count>=mini_clust_size]
listmatch=vector(mode="list")
for(clust in clust_filter){
  ind_clust=which(nodetab[,"__mclCluster"]==clust)
  nodenames=nodetab[ind_clust,"name"]
  ppmvec=str_extract_all(string=nodenames,pattern="\\d\\.\\d+$")
  # search on gissmo
  apipara[["rs"]]=paste0(ppmvec,collapse=" ")
  requests<-import("requests")
  res=requests$get(peaksearch_url,params=apipara)
  reslist=res$json()
  lenmatch=sapply(reslist,function(x){
    length(x[["Val"]])
  })
  mathcind_ord=order(lenmatch,decreasing=TRUE)
  ratio_ind=mathcind_ord[lenmatch[mathcind_ord]/length(ppmvec)>=mini_match_prop | lenmatch[mathcind_ord]>mini_match_num]#the match will be visualized if one of the condition is true
  order_ind=mathcind_ord[1:maxmatch]
  peakind=intersect(ratio_ind,order_ind)

  cleanlist=sapply(reslist[peakind],simplify=FALSE,function(x){
    x[1:4]
  })
  matched_peak_tab=as.data.frame(Reduce(rbind.data.frame,cleanlist))
  # visualize in MATLAB
  for(rowi in 1:nrow(matched_peak_tab)){
    message(matched_peak_tab[rowi,"Entry_ID"])
    script=paste0('matchnode={\'',paste0(nodenames,collapse="\',\'"),'\'}; ',
                 'id=\'',matched_peak_tab[rowi,"Entry_ID"],'\'; ',
                 'sampi=',show_samp_i,'; ',
                 'paralist.datadir=\'',pardir,'data/unshared/\'; ',
                 'paralist.gissmodir=\'',pardir,'data/shared/\'; ',
                 'paralist.match_dir=\'',datadir,'\'; ',
                 'flag_run=spec_vis_gissmo(id,sampi,matchnode,paralist); ',
                 'if flag_run~=0; exit; end; ',
                 'fig=gca; saveas(fig,\'anno_show_',clust,'_',rowi,'.fig\'); close all; ',
                 'exit')
    cat(script,file="temp_script.m")
    system('matlab -r \"run(\'./temp_script.m\'); exit\"');
  }
  listmatch[[clust]]=list(nodes=nodenames,matched=matched_peak_tab)
}
save(listmatch,file="stored_matched_cluster.RData")
