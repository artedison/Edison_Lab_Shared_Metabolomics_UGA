rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
dir=""
shelltempt=paste0(dir,"submit.sh")#submit.sh need to be a proper shell script for the HPC
nnetworks=60
for(starti in seq(nnetworks)){
  runshell=paste0("runshell",starti,".sh")
  system(paste0("cp ",shelltempt," ",runshell))
  lines=readLines(runshell)
  chline_ind=str_which(string=lines,pattern="^time")
  lineend=str_extract_all(string=lines[chline_ind],pattern="\\w+>>.*$")[[1]]
  lines[chline_ind]=paste("time Rscript simu_network_testcluster.R",starti,lineend,sep=" ")
  cat(lines,file=runshell,sep="\n")
  submitcommand=paste0("sbatch ",runshell)
  print(submitcommand)
  system(submitcommand)
}
