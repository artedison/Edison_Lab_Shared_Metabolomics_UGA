# The script to load bootstraping constructed network and highlight consistent edges in the original network
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
# the user will need to modify the path here for local run
comp="/Users/yuewu/"
pardir=paste0(comp,"Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/");
# OR
# pardir=paste0(getwd(),"/")
datadir=paste0(pardir,"result_data/")
resdir=paste0(pardir,"result/clust_bootstrapping/")
setwd(resdir)
load(paste0(datadir,"connection_tab.RData"))
cytoscapePing()
cytoscapeVersionInfo()
file.copy(paste0(datadir,"net.causalkinetix.full.pairwise.cys"),paste0(resdir,"net.causalkinetix.full.pairwise.copy.cys"))
openSession(paste0(resdir,"net.causalkinetix.full.pairwise.copy.cys"))
setCurrentView("net.CausalKinetiXpeakwise.txt--clustered")
edgetab=getTableColumns(table="edge")
edgesize=dim(edgetab)
edgetab$"bootstrap_flag"=rep("No_consistent",times=edgesize[1])
filttab=listtab[["res_peakwise"]]
nedge=dim(filttab)[1]
nodenames=colnames(filttab)
thres_bootstrap=40#0.4
control_points=c('No_consistent','Consistent')
colors=c('#808080','#FF0000')
style="BioPAX_0_0"
for(edgei in seq(nedge-1)){
  for(edgej in seq(from=edgei+1,by=1,to=nedge)){
    if(filttab[edgei,edgej]>=thres_bootstrap){
      edgename=c(paste0(nodenames[edgei]," (interacts with) ",nodenames[edgej]),paste0(nodenames[edgej]," (interacts with) ",nodenames[edgei]))
      matchededge=c(str_which(string=edgetab[,"name"],pattern=fixed(edgename[1])),str_which(string=edgetab[,"name"],pattern=fixed(edgename[2])))
      if(length(matchededge)>1){
        stop("duplicated edge")
      }
      edgetab[matchededge,"bootstrap_flag"]="Consistent"
    }
  }
}
loadTableData(edgetab,table="edge",data.key.column="shared name",table.key.column="shared name")
setVisualStyle(style)
setEdgeColorMapping("bootstrap_flag",control_points,colors,mapping.type="discrete",style.name=style)
setEdgeLineWidthMapping("bootstrap_flag",control_points,c(1.0,5.0),mapping.type="discrete",style.name=style)
setEdgeOpacityMapping("bootstrap_flag",control_points,c(50,255),mapping.type="discrete",style.name=style)
saveSession("connection_tab.highlighted.cys")
