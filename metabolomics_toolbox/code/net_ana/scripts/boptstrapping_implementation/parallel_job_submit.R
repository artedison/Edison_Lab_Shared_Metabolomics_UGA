# batch submision script for CausalKinetiX.cluster.bootstrapping.R
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
dir=""
shelltempt=paste0(dir,"submit.sh")#submit.sh need to be a proper shell script for the HPC
nbootstrp=100
nbootstrapstep=10
for(starti in seq(from=1,by=nbootstrapstep,to=nbootstrp)){
  bootstrapseq=c(starti,starti+nbootstrapstep-1)
  runshell=paste0("runshell",starti,".sh")
  system(paste0("cp ",shelltempt," ",runshell))
  lines=readLines(runshell)
  chline_ind=str_which(string=lines,pattern="^time")
  lineend=str_extract_all(string=lines[chline_ind],pattern="\\w+>>.*$")[[1]]
  lines[chline_ind]=paste("time Rscript CausalKinetiX.cluster.bootstrapping.R",bootstrapseq[1],bootstrapseq[2],lineend,sep=" ")
  cat(lines,file=runshell,sep="\n")
  submitcommand=paste0("qsub ",runshell)
  print(submitcommand)
  system(submitcommand)
}
