## this script is for evaluate and plot the simulated time series spectra
## 1. intensity of ridges for several example compound
## 2. linearity/correlation of intensity and ppm
## 3. affect on noise level on ridge tracing
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(R.matlab)
require(ggplot2)
require(xml2)
require(reshape2)
require(dplyr)
require(scales)
require(plotly)
require(DescTools)
require(cowplot)
## load real concentration
comp="/Users/yuewu/"
dir=paste0(comp,"/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/result/")
dir.data=paste0(comp,"/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/data/")
dir.match=paste0(dir,"peakmatching_simulated/")
dir.r.res=paste0(dir,"peakmatching_simulated/R_result/")
matres=readMat(paste0(dir.match,"eval.ridtracing.mat"))
rawdata=readMat(paste0(dir,"simulateddata/simulated.timeseries.ph.complex.more.addunknown.mat"))
peaklistraw=readMat(paste0(dir.data,"spectral.real.list.withshift.mat"))
ppmrealtab=peaklistraw[["strdatalist"]]## the table for ppm list and conerter of peak intensity
comptruevallist=rawdata[["truevallist"]]## list of true intensity
compds=c("Acetate","Formate","ethanol","Glycerol","Choline","alanine","Uridine","Leucine","Valine","DSS","Butanol","Caffeine","Serine","Purine","unknown1","unknown2","unknown3","unknown4")
ncompds=length(compds);
datasources="manual"#c("manual","automatic")
inforused=c("intensity","ppm","time")
samples=4
ncol=3##number of colume in combined plot
ppmthred=0.00001
shiftcompd=c("Acetate","Formate","unknown1","unknown2","unknown3","unknown4")
samplelevels=seq(samples)
dataset=matres[["list.match"]]##list of matched peaks
## convert to a dataframe
datasource.indata=attr(dataset,'dimnames')[[1]]## this part is to make sure the using order for data from matlab is correct
##this is a table for peak: compound name, ridge id, datasources, diferent simulated sample, ppm, intensity, time, index for different ridge in same compound
datalist=list("compd"=c(),"ridid"=c(),"datasource"=c(),"samplelevel"=c(),
              "ppm"=c(),"intensity"=c(),"time"=c(),"ridgeindex"=c(),"matchppm"=c())
##reformat data reform the matched peak into a dataframe
for(datasource in datasources){## for each data source
  datasource.ind=which(datasource==datasource.indata)
  datasamp=dataset[[datasource.ind]]
  # datasamp.indata=attr(datasamp,'dimnames')[[1]]
  for(samplelevel in samplelevels){##different sample for different complex level
    # datasamp.ind=which(samplelevel==datasamp.indata)
    # datacompd=datasamp[[datasamp.ind]]
    datacompd=datasamp[[samplelevel]][[1]]
    datacompd.indata=attr(datacompd,'dimnames')[[1]]
    for(compd in datacompd.indata){## for each compound
      datacompd.ind=which(datacompd.indata==compd)
      datainside=datacompd[[datacompd.ind]]
      for(ridi in seq(length(datainside))){##for each ridge
        datarid=datainside[[ridi]][[1]]
        datarid.indata=attr(datarid,'dimnames')[[1]]
        datares.ind=which(datarid.indata=="result")
        datares=datarid[[datares.ind]]
        matchppm=datarid[[which(datarid.indata=="match")]][1]
        datares.indata=attr(datares,'dimnames')[[1]]
        lengthrid=0;
        for(field in inforused){## for each struct field
          field.ind=which(field==datares.indata)
          datalist[[field]]=c(datalist[[field]],datares[[field.ind]])
          lengthrid=length(datares[[field.ind]]);
        }
        ppmmean=mean(datares[[which(datares.indata=="ppm")]])
        ppmmean=formatC(ppmmean,format="f",digits=3)
        datalist[["ridgeindex"]]=c(datalist[["ridgeindex"]],rep(ppmmean,times=lengthrid))
        ridgeid=str_replace_all(string=paste(datasource,samplelevel,compd,ridi,sep="_"),pattern="\\s",replacement="_")
        datalist[["ridid"]]=c(datalist[["ridid"]],rep(ridgeid,times=lengthrid))
        for(otherfield in setdiff(names(datalist),c(inforused,"ridid","ridgeindex"))){
          datalist[[otherfield]]=c(datalist[[otherfield]],rep(get(otherfield),times=lengthrid))
        }
      }
    }
  }
}
riddataframe=as.data.frame(datalist)
colnames(riddataframe)=names(datalist)
###numeric conversion
for(col in c("ppm","intensity","time")){
  riddataframe[,col]=as.numeric(riddataframe[,col])
}
###factor conversion
for(col in c("ridid")){
  riddataframe[,col]=as.factor(riddataframe[,col])
}
save(riddataframe,file=paste0(dir.r.res,"eval.ridtracing.RData"))
##scale the intensity for each ridges
load(paste0(dir.r.res,"eval.ridtracing.RData"))
##add true value
# strdataannotname=c("Acetate","Formate","ethanol","Glycerol","Choline","alanine","Uridine","Leucine","Valine")
riddataframe$trueval=rep(NA,times=dim(riddataframe)[1])
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    comptruevallisttemp=comptruevallist[[ceiling(samplelevelhere/2)]]
    compdslist=attr(comptruevallisttemp[[1]],'dimnames')[[1]]
    for(compdhere in compdslist){
      numvec=comptruevallisttemp[[1]][compdhere,,][[1]][1,]
      ind=(riddataframe[,"datasource"]==datasourcehere & riddataframe[,"samplelevel"]==as.character(samplelevelhere) & riddataframe[,"compd"]==compdhere)
      maxconc=max(numvec)
      # tab=riddataframe[ind,]
      if(compdhere=="DSS"){
        numvec=rep(1,times=length(numvec))
      }else{
        numvec=rescale(numvec)
      }
      factortable=ppmrealtab[,,1][[compdhere]]
      for(ridhere in unique(riddataframe[ind,"ridid"])){
        indrid=(riddataframe[,"ridid"]==ridhere)
        timevec=riddataframe[riddataframe[,"ridid"]==ridhere,"time"]
        riddataframe[indrid,"trueval"]=numvec[timevec]##time works as a index to make the true value and measured value corresponding to each other
        peakppm=unique(riddataframe[indrid,"matchppm"])
        if(!is.null(factortable)){
          matind=abs(factortable[,1]-peakppm)<ppmthred
          if(length(which(matind))>1){
            stop("overlapped ppm")
          }
          range=c(0,factortable[matind,2]*maxconc)
          if(compdhere%in%shiftcompd){
            range=c(0,1*maxconc)
          }
          scaledintensity=rescale(riddataframe[riddataframe[,"ridid"]==ridhere,"intensity"],from=range,to=c(0,1))
        }else{
          scaledintensity=rescale(riddataframe[riddataframe[,"ridid"]==ridhere,"intensity"])
        }
        riddataframe[indrid,"intensity"]=scaledintensity
      }
    }
  }
}
riddataframe[,"matchppm"]=as.character(riddataframe[,"matchppm"])
##mutiple ridges for same compound
nochangeconc=c("DSS","unknown1","unknown2","unknown3")##compound with no changes in quantity
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    plist=vector(mode="list")
    for(compdhere in compds){
      tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")) & compd==get("compdhere"))
      if(dim(tab)[1]==0){
        next
      }
      name=paste0(compdhere,"_samplelevel",samplelevelhere,"_",datasourcehere)
      # p<-ggplot(data=tab,aes(time,intensity,color=ridid))+
      #       geom_point(alpha=0.75)+
      #       xlab("time(h)")+
      #       ylab("scaled_intensity")+
      #       ggtitle(name)+
      #       guides(color=guide_legend(title="ridge_ppm"))
      #       theme_bw()+
      #       theme(plot.title=element_text(hjust=0.5))
      # ggsave(plot=p,file=paste0(dir.r.res,name,".timetraj.pdf"))
      p<-ggplot(data=tab,aes(trueval,intensity,color=matchppm))+
            geom_point(alpha=0.75)+
            xlab("true value")+
            ylab("scaled_intensity")+
            ggtitle(name)+
            guides(color=guide_legend(title="ridge_ppm"))
            theme_bw()+
            theme(plot.title=element_text(hjust=0.5))
      ggsave(plot=p,file=paste0(dir.r.res,name,".trueval_esitmatedval_traj.pdf"))
      if(!compdhere%in%nochangeconc){
        p=p+theme(legend.position="none",plot.title=element_blank(),
                  axis.title.x=element_blank(),axis.title.y=element_blank(),
                  axis.text.x=element_text(angle=90,hjust=1))
        plist[[compdhere]]=p
      }
    }
    ##plot a combined figure
    p<-plot_grid(plotlist=plist,labels=c(names(plist)),
            label_size=10,ncol=ncol,label_x=0.3)
    ggsave(plot=p,filename=paste0(dir.r.res,"samplelevel",samplelevelhere,"_",datasourcehere,"trueval_esitmatedval_traj_sum.pdf"),
            limitsize=FALSE)
  }
}
##add true ppm value
ppmchtable=data.frame(compd=c("Acetate","Acetate","Formate","Formate","unknown2","unknown3","unknown4"),complex=c(1,2,1,2,2,2,2))
riddataframe$trueppm=rep(0,times=dim(riddataframe)[1])
for(datasourcehere in datasources){
  for(i in seq(dim(ppmchtable)[1])){
    setind=ppmchtable[i,]
    compdhere=setind[1,1]
    dimnames_compd=attr(ppmrealtab,'dimnames')[[1]]
    ppmreal=c(ppmrealtab[[which(dimnames_compd==compdhere)]])
    complexlevelhere=setind[1,2]
    samplelevelheres=c(complexlevelhere*2,complexlevelhere*2-1)
    for(samplelevelhere in samplelevelheres){
      indrid=(riddataframe[,"compd"]==compdhere&riddataframe[,"samplelevel"]==as.character(samplelevelhere)&riddataframe[,"datasource"]==datasourcehere)
      timevec=riddataframe[indrid,"time"]
      riddataframe[indrid,"trueppm"]=ppmreal[timevec]##time works as a index to make the true value and measured value corresponding to each other
    }
  }
}
save(riddataframe,file=paste0(dir.r.res,"eval.ridtracing.withtrue.RData"))
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    plist=vector(mode="list")
    for(compdhere in compds){
      tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")) & compd==get("compdhere"))
      if(dim(tab)[1]==0){
        next
      }
      name=paste0(compdhere,"_samplelevel",samplelevelhere,"_",datasourcehere)
      p<-ggplot(data=tab,aes(trueppm,ppm,color=matchppm))+
            geom_point(alpha=0.75)+
            xlab("true ppm")+
            ylab("traced ppm")+
            ggtitle(name)+
            guides(color=guide_legend(title="ridge_ppm"))
            theme_bw()+
            theme(plot.title=element_text(hjust=0.5))
      if(compdhere=="Acetate"){
        p=p+scale_x_continuous(limits=c(1.9,2.05))+scale_y_continuous(limits=c(1.9,2.05))
      }
      ggsave(plot=p,file=paste0(dir.r.res,name,".trueppm_esitmatedppm_traj.pdf"))
      if(compdhere %in% ppmchtable[,"compd"]){
        p=p+theme(legend.position="none",plot.title=element_blank(),
                  axis.title.x=element_blank(),axis.title.y=element_blank(),
                  axis.text.x=element_text(angle=90,hjust=1))
        plist[[compdhere]]=p
      }
    }
    ##plot a combined figure
    p<-plot_grid(plotlist=plist,labels=c(names(plist)),
            label_size=10,ncol=ncol,label_x=0.3)
    ggsave(plot=p,filename=paste0(dir.r.res,"samplelevel",samplelevelhere,"_",datasourcehere,"trueppm_esitmatedppm_traj_sum.pdf"),
            limitsize=FALSE)
  }
}

###global plot for absolute unscaled value for ppm and intensity for all simulated dataset
####preprocessing
load(paste0(dir.r.res,"eval.ridtracing.RData"))
riddataframe$trueval=rep(NA,times=dim(riddataframe)[1])
riddataframe$trueppm=rep(NA,times=dim(riddataframe)[1])
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    comptruevallisttemp=comptruevallist[[ceiling(samplelevelhere/2)]]
    compdslist=attr(comptruevallisttemp[[1]],'dimnames')[[1]]
    for(compdhere in compdslist){
      numvec=comptruevallisttemp[[1]][compdhere,,][[1]][1,]
      ind=(riddataframe[,"datasource"]==datasourcehere & riddataframe[,"samplelevel"]==as.character(samplelevelhere) & riddataframe[,"compd"]==compdhere)
      maxconc=max(numvec)
      # tab=riddataframe[ind,]
      if(compdhere=="DSS"){
        numvec=rep(1,times=length(numvec))
      }else{
        numvec=rescale(numvec)
      }
      factortable=ppmrealtab[,,1][[compdhere]]
      for(ridhere in unique(riddataframe[ind,"ridid"])){
        indrid=(riddataframe[,"ridid"]==ridhere)
        timevec=riddataframe[riddataframe[,"ridid"]==ridhere,"time"]
        riddataframe[indrid,"trueval"]=numvec[timevec]*maxconc##time works as a index to make the true value and measured value corresponding to each other
        peakppm=unique(riddataframe[indrid,"matchppm"])
        if(!is.null(factortable)){
          matind=abs(factortable[,1]-peakppm)<ppmthred
          if(length(which(matind))>1){
            stop("overlapped ppm")
          }
          range=c(0,factortable[matind,2]*maxconc)
          if(compdhere%in%shiftcompd){
            range=c(0,1*maxconc)
          }
          scaledintensity=rescale(riddataframe[riddataframe[,"ridid"]==ridhere,"intensity"],from=range,to=c(0,1))
          ##ppm that not shift
          if(!compdhere%in%shiftcompd){
            riddataframe[indrid,"trueppm"]=factortable[matind,1]
          }
          if(compdhere=="unknown1"){
            riddataframe[indrid,"trueppm"]=factortable[1,1]
          }
        }else{
          scaledintensity=rescale(riddataframe[riddataframe[,"ridid"]==ridhere,"intensity"])
        }
        riddataframe[indrid,"intensity"]=scaledintensity*maxconc
      }
    }
  }
}
riddataframe[,"matchppm"]=as.character(riddataframe[,"matchppm"])
ppmchtable=data.frame(compd=c("Acetate","Acetate","Formate","Formate","unknown2","unknown3","unknown4"),complex=c(1,2,1,2,2,2,2))
for(datasourcehere in datasources){
  for(i in seq(dim(ppmchtable)[1])){
    setind=ppmchtable[i,]
    compdhere=setind[1,1]
    dimnames_compd=attr(ppmrealtab,'dimnames')[[1]]
    ppmreal=c(ppmrealtab[[which(dimnames_compd==compdhere)]])
    complexlevelhere=setind[1,2]
    samplelevelheres=c(complexlevelhere*2,complexlevelhere*2-1)
    for(samplelevelhere in samplelevelheres){
      indrid=(riddataframe[,"compd"]==compdhere&riddataframe[,"samplelevel"]==as.character(samplelevelhere)&riddataframe[,"datasource"]==datasourcehere)
      timevec=riddataframe[indrid,"time"]
      riddataframe[indrid,"trueppm"]=ppmreal[timevec]##time works as a index to make the true value and measured value corresponding to each other
    }
  }
}
save(riddataframe,file=paste0(dir.r.res,"eval.ridtracing.global.RData"))
#### for intensity
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")))
    name=paste0("global_intensity_samplelevel",samplelevelhere,"_",datasourcehere)
    color=factor(ifelse(tab[,"compd"] %in% c("ethanol"),"ethanol","others"))
    p<-ggplot(data=tab,aes(trueval,intensity,color=color))+#
          geom_point(alpha=0.1)+
          scale_colour_manual(values=c("ethanol"="red","others"="grey"))+
          xlab("true intensity")+
          ylab("traced intensity")+
          xlim(-20,3200)+
          ylim(-20,3200)+
          ggtitle(name)+
          # guides(color=guide_legend(title="ridge_ppm"))+
          theme_bw()+
          theme(plot.title=element_text(hjust=0.5))+
          theme(legend.position="none",
                axis.text.x=element_text(angle=90,hjust=1))
    ggsave(plot=p,filename=paste0(dir.r.res,"samplelevel",samplelevelhere,"_",datasourcehere,"trueintensity_esitmatedintensity_traj_global_ethanol.pdf"),
           limitsize=FALSE)
  }
}
#### for ppm
shiftcopd=unique(ppmchtable[,1])
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")))
    name=paste0("global_ppm_samplelevel",samplelevelhere,"_",datasourcehere)
    color=factor(ifelse(tab[,"compd"] %in% shiftcopd,"pH_shift","pH_stable"))
    p<-ggplot(data=tab,aes(trueppm,ppm,color=color))+
          geom_point(alpha=0.1)+
          scale_colour_manual(values=c("pH_shift"="red","pH_stable"="grey"))+
          xlab("true ppm")+
          ylab("traced ppm")+
          xlim(-0.5,10)+
          ylim(-0.5,10)+
          ggtitle(name)+
          # guides(color=guide_legend(title="ridge_ppm"))+
          theme_bw()+
          theme(plot.title=element_text(hjust=0.5))+
          theme(legend.position="none",
                axis.text.x=element_text(angle=90,hjust=1))
    ggsave(plot=p,filename=paste0(dir.r.res,"samplelevel",samplelevelhere,"_",datasourcehere,"trueppm_esitmatedppm_traj_global.pdf"),
           limitsize=FALSE)
  }
}

###global plot for absolute unscaled value for ppm and intensity for all simulated dataset
## this block use dss peak to normalize to calculate estimated compound concentration
####preprocessing
load(paste0(dir.r.res,"eval.ridtracing.RData"))
riddataframe$trueval=rep(NA,times=dim(riddataframe)[1])
riddataframe$trueppm=rep(NA,times=dim(riddataframe)[1])
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    comptruevallisttemp=comptruevallist[[ceiling(samplelevelhere/2)]]
    compdslist=attr(comptruevallisttemp[[1]],'dimnames')[[1]]
    ##reference information
    refcompd="DSS"
    normal.struc=vector(mode="list")##DSS
    numvec=comptruevallisttemp[[1]][refcompd,,][[1]][1,]
    ind=(riddataframe[,"datasource"]==datasourcehere & riddataframe[,"samplelevel"]==as.character(samplelevelhere) & riddataframe[,"compd"]==refcompd)
    refconc=numvec[1]
    factorref=ppmrealtab[,,1][[refcompd]][1,]##information for the reference peak of DSS
    normal.struc[["realconcentr_norm"]]=refconc*factorref[2]
    peakppmvec=unique(riddataframe[ind,"matchppm"])
    matindref=abs(peakppmvec-factorref[1])<ppmthred
    indridref=(ind&riddataframe[,"matchppm"]==peakppmvec[matindref])
    normal.struc[["timevec"]]=riddataframe[indridref,"time"]
    normal.struc[["peakvec"]]=riddataframe[indridref,"intensity"]
    for(compdhere in compdslist){
      numvec=comptruevallisttemp[[1]][compdhere,,][[1]][1,]
      ind=(riddataframe[,"datasource"]==datasourcehere & riddataframe[,"samplelevel"]==as.character(samplelevelhere) & riddataframe[,"compd"]==compdhere)
      maxconc=max(numvec)
      # tab=riddataframe[ind,]
      # if(compdhere=="DSS"){
      #   numvec=rep(1,times=length(numvec))
      # }else{
      #   numvec=rescale(numvec)
      # }
      factortable=ppmrealtab[,,1][[compdhere]]
      for(ridhere in unique(riddataframe[ind,"ridid"])){
        indrid=(riddataframe[,"ridid"]==ridhere)
        timevec=riddataframe[riddataframe[,"ridid"]==ridhere,"time"]
        riddataframe[indrid,"trueval"]=numvec[timevec]##time works as a index to make the true value and measured value corresponding to each other
        peakppm=unique(riddataframe[indrid,"matchppm"])
        if(!is.null(factortable)){
          matind=abs(factortable[,1]-peakppm)<ppmthred
          if(length(which(matind))>1){
            stop("overlapped ppm")
          }
          # range=c(0,factortable[matind,2]*maxconc)
          # if(compdhere%in%shiftcompd){
          #   range=c(0,1*maxconc)
          # }
          # scaledintensity=rescale(riddataframe[riddataframe[,"ridid"]==ridhere,"intensity"],from=range,to=c(0,1))
          ridscal=factortable[matind,2]
          if(compdhere%in%shiftcompd){
            ridscal=1
          }
          scaledintensity=riddataframe[indrid,"intensity"]*normal.struc[["realconcentr_norm"]]/normal.struc[["peakvec"]][timevec]/ridscal
          ##ppm that not shift
          if(!compdhere%in%shiftcompd){
            riddataframe[indrid,"trueppm"]=factortable[matind,1]
          }
          if(compdhere=="unknown1"){
            riddataframe[indrid,"trueppm"]=factortable[1,1]
          }
        }else{
          stop("problem with peak evaluation")
          # scaledintensity=rescale(riddataframe[riddataframe[,"ridid"]==ridhere,"intensity"])
        }
        riddataframe[indrid,"intensity"]=scaledintensity#*maxconc
      }
    }
  }
}
riddataframe[,"matchppm"]=as.character(riddataframe[,"matchppm"])
ppmchtable=data.frame(compd=c("Acetate","Acetate","Formate","Formate","unknown2","unknown3","unknown4"),complex=c(1,2,1,2,2,2,2))
for(datasourcehere in datasources){
  for(i in seq(dim(ppmchtable)[1])){
    setind=ppmchtable[i,]
    compdhere=setind[1,1]
    dimnames_compd=attr(ppmrealtab,'dimnames')[[1]]
    ppmreal=c(ppmrealtab[[which(dimnames_compd==compdhere)]])
    complexlevelhere=setind[1,2]
    samplelevelheres=c(complexlevelhere*2,complexlevelhere*2-1)
    for(samplelevelhere in samplelevelheres){
      indrid=(riddataframe[,"compd"]==compdhere&riddataframe[,"samplelevel"]==as.character(samplelevelhere)&riddataframe[,"datasource"]==datasourcehere)
      timevec=riddataframe[indrid,"time"]
      riddataframe[indrid,"trueppm"]=ppmreal[timevec]##time works as a index to make the true value and measured value corresponding to each other
    }
  }
}
save(riddataframe,file=paste0(dir.r.res,"eval.ridtracing.global.dssref.RData"))
#### for intensity
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")))
    name=paste0("global_intensity_samplelevel",samplelevelhere,"_",datasourcehere)
    # color=factor(ifelse(tab[,"compd"] %in% c("ethanol"),"ethanol","others"))
    # color=factor(tab[,"compd"])
    p<-ggplot(data=tab,aes(trueval,intensity))+#,color=color
          geom_point(alpha=0.1)+
          # scale_colour_manual(values=c("ethanol"="red","others"="grey"))+
          xlab("true intensity")+
          ylab("traced intensity")+
          xlim(-100,4400)+
          ylim(-100,4400)+
          ggtitle(name)+
          # guides(color=guide_legend(title="ridge_ppm"))+
          theme_bw()+
          theme(plot.title=element_text(hjust=0.5))+
          theme(legend.position="none",
                axis.text.x=element_text(angle=90,hjust=1))
    # p<-ggplotly(p)
    ggsave(plot=p,filename=paste0(dir.r.res,"samplelevel",samplelevelhere,"_",datasourcehere,"trueintensity_esitmatedintensity_traj_global.dssref.pdf"),
           limitsize=FALSE)
  }
}
#### for ppm
shiftcopd=unique(ppmchtable[,1])
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")))
    name=paste0("global_ppm_samplelevel",samplelevelhere,"_",datasourcehere)
    color=factor(ifelse(tab[,"compd"] %in% shiftcopd,"pH_shift","pH_stable"))
    p<-ggplot(data=tab,aes(trueppm,ppm,color=color))+
          geom_point(alpha=0.1)+
          scale_colour_manual(values=c("pH_shift"="red","pH_stable"="grey"))+
          xlab("true ppm")+
          ylab("traced ppm")+
          xlim(-0.5,10)+
          ylim(-0.5,10)+
          ggtitle(name)+
          # guides(color=guide_legend(title="ridge_ppm"))+
          theme_bw()+
          theme(plot.title=element_text(hjust=0.5))+
          theme(legend.position="none",
                axis.text.x=element_text(angle=90,hjust=1))
    ggsave(plot=p,filename=paste0(dir.r.res,"samplelevel",samplelevelhere,"_",datasourcehere,"trueppm_esitmatedppm_traj_globaldssref.pdf"),
           limitsize=FALSE)
  }
}
##highlight non-overlap region
# ## good region for quantificaiton region selected only based on sample 4
# regionspeaks=struct();
## regionspeaks.Acetate=[the only peak overlap with others]
## regionspeaks.Formate=[8.4 8.5];
## regionspeaks.ethanol=[1.15 1.185];%3%%
## regionspeaks.Choline=[3.18 3.22];
## regionspeaks.alanine=[1.45 1.48];%%
# regionspeaks.Uridine=[3.89 3.9],[4.1 4.14] [4.3 4.4] [5.86 5.92]
## regionspeaks.Leucine=[0.944 0.954];%%
## regionspeaks.Valine=[1.02 1.04] [2.2 2.4];%2%%
## regionspeaks.DSS=[-0.01 0.01];%%
## Glycerol [3.5 3.6]
# Butanol [0.85 0.92] [1.3 1.4] [3.57 3.61]
# Caffeine [all peaks overalp with others]
# Serine [3.93 3.99]
# Purine [8.5 9.1]
# unkown1 [the only peak overlap with others]
# unkown2 [4.19 4.26]
# unkown3 [4.84 4.91]
# unkown4 [3.34 3.46]
region.list=vector(mode="list")
region.list[["Formate"]]=c(8.4,8.5)
region.list[["ethanol"]]=c(1.15,1.185)
region.list[["Choline"]]=c(3.18,3.22)
region.list[["alanine"]]=c(1.45,1.48)
region.list[["Uridine1"]]=c(3.89,3.9)
region.list[["Uridine2"]]=c(4.1,4.14)
region.list[["Uridine3"]]=c(4.3,4.4)
region.list[["Uridine4"]]=c(5.86,5.92)
region.list[["Leucine"]]=c(0.944,0.954)
region.list[["Valine1"]]=c(1.02,1.04)
region.list[["Valine2"]]=c(2.2,2.4)
region.list[["DSS"]]=c(-0.01,0.01)
region.list[["Glycerol"]]=c(3.5,3.6)
region.list[["Butanol1"]]=c(0.85,0.92)
region.list[["Butanol2"]]=c(1.3,1.4)
region.list[["Butanol3"]]=c(3.57,3.61)
region.list[["Serine"]]=c(3.93,3.99)
region.list[["Purine"]]=c(8.5,9.1)
region.list[["unkown2"]]=c(4.19,4.26)
region.list[["unkown3"]]=c(4.84,4.91)
region.list[["unkown4"]]=c(3.34,3.46)
hightlightcompds=names(region.list)
# listres=combn(names(region.list),2,simplify=FALSE,function(x){
#           name1=x[[1]]
#           name2=x[[2]]
#           c(name1,name2,x[[1]] %overlaps% x[[2]])
#         }) #
# tabres=Reduce("rbind",listres)
# combn(names(region.list),2,function(x){
#           x[[1]] %overlaps% x[[2]]
#         }) #
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")))
    name=paste0("global_intensity_samplelevel",samplelevelhere,"_",datasourcehere)
    color=rep("overlap",times=dim(tab)[1])
    for(hightlightcompd in hightlightcompds){
      range=region.list[[hightlightcompd]]
      hightlightcompd=str_replace_all(string=hightlightcompd,pattern="\\d+",replacement="")
      matchedppm=as.numeric(tab[,"matchppm"])
      indmatch=(tab[,"compd"]==hightlightcompd&matchedppm>range[1]&matchedppm<range[2])
      color[indmatch]="nonoverlap"
    }
    color=factor(color)
    shape=factor(tab[,"compd"])
    p<-ggplot(data=tab,aes(trueval,intensity,color=color,shape=shape))+#
          geom_point(alpha=0.7,size=0.3)+
          # geom_abline(slope=1,intercept=0,color="black")+
          geom_segment(aes(x=0,xend=2000,y=0,yend=2000),color="black")+
          scale_colour_manual(values=c("nonoverlap"="red","overlap"="grey"))+
          xlab("true intensity")+
          ylab("traced intensity")+
          xlim(-100,4400)+
          ylim(-100,4400)+
          ggtitle(name)+
          # guides(color=guide_legend(title="ridge_ppm"))+
          theme_bw()+
          theme(plot.title=element_text(hjust=0.5))+
          theme(legend.position="none",
                axis.text.x=element_text(angle=90,hjust=1))
    # p<-ggplotly(p)
    # p
    ggsave(plot=p,filename=paste0(dir.r.res,"samplelevel",samplelevelhere,"_",datasourcehere,"trueintensity_esitmatedintensity_traj_global.dssref.nonoverlaphightlight.pdf"),
           limitsize=FALSE)
  }
}
##hightlight non-overlap and also color specific peaks
##unique(riddataframe[riddataframe[,"matchppm"]=="ppm","ridid"])##check
Peakofinterest=c("alanine_1"="3.76133","Serine_2"="3.91635","Glycerol_11"="3.77877","ethanol_4"="3.62311")
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")))
    name=paste0("global_intensity_samplelevel",samplelevelhere,"_",datasourcehere)
    color=rep("overlap",times=dim(tab)[1])
    for(hightlightcompd in hightlightcompds){
      range=region.list[[hightlightcompd]]
      hightlightcompd=str_replace_all(string=hightlightcompd,pattern="\\d+",replacement="")
      matchedppm=as.numeric(tab[,"matchppm"])
      indmatch=(tab[,"compd"]==hightlightcompd&matchedppm>range[1]&matchedppm<range[2])
      color[indmatch]="nonoverlap"
    }
    for(peakname in names(Peakofinterest)){
      peak=Peakofinterest[[peakname]]
      cmp_color_ind=(tab[,"matchppm"]==peak)
      color[cmp_color_ind]=peakname
    }
    color=factor(color)
    # shape=factor(tab[,"compd"])
    p<-ggplot(data=tab,aes(trueval,intensity,color=color))+#,shape=shape
          geom_point(alpha=0.8,size=1.2)+
          # geom_abline(slope=1,intercept=0,color="black")+
          geom_segment(aes(x=0,xend=2000,y=0,yend=2000),color="black")+
          scale_colour_manual(values=c("nonoverlap"="red","overlap"="grey","alanine_1"="blue","Serine_2"="darkturquoise","Glycerol_11"="goldenrod","ethanol_4"="darkviolet"))+
          xlab("true concentration")+
          ylab("normlaized estimated concentration")+
          xlim(-100,4400)+
          ylim(-100,4400)+
          ggtitle(name)+
          guides(color=guide_legend(title="ridge_ppm"))+
          # theme_bw()+
          theme(plot.title=element_text(hjust=0.5))+
          theme(legend.position="right",plot.title=element_blank(),
                axis.title.x=element_blank(),axis.title.y=element_blank(),
                axis.text.x=element_text(angle=90,hjust=1))
    # p<-ggplotly(p)
    # p
    ggsave(plot=p,filename=paste0(dir.r.res,"samplelevel",samplelevelhere,"_",datasourcehere,"trueintensity_esitmatedintensity_traj_global.dssref.nonoverlaphightlight.addspec.pdf"),
           limitsize=FALSE)
  }
}

##count ridges
length(unique(tab[color=="nonoverlap","matchppm"]))
length(unique(tab[,"matchppm"]))

##plotting individual intensity
nochangeconc=c("DSS","unknown1","unknown2","unknown3")##compound with no changes in quantity
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    plist=vector(mode="list")
    for(compdhere in compds){
      tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")) & compd==get("compdhere"))
      if(dim(tab)[1]==0){
        next
      }
      name=paste0(compdhere,"_samplelevel",samplelevelhere,"_",datasourcehere)
      p<-ggplot(data=tab,aes(trueval,intensity,color=matchppm))+
            geom_point(alpha=0.75)+
            xlab("true value")+
            ylab("scaled_intensity")+
            ggtitle(name)+
            guides(color=guide_legend(title="ridge_ppm"))
            theme_bw()+
            theme(plot.title=element_text(hjust=0.5))
      ggsave(plot=p,file=paste0(dir.r.res,name,".trueval_esitmatedval_traj_dssref.pdf"))
      if(!compdhere%in%nochangeconc){
        p=p+theme(legend.position="none",plot.title=element_blank(),
                  axis.title.x=element_blank(),axis.title.y=element_blank(),
                  axis.text.x=element_text(angle=90,hjust=1))
        plist[[compdhere]]=p
      }
    }
    ##plot a combined figure
    p<-plot_grid(plotlist=plist,labels=c(names(plist)),
            label_size=10,ncol=ncol,label_x=0.3)
    ggsave(plot=p,filename=paste0(dir.r.res,"samplelevel",samplelevelhere,"_",datasourcehere,"trueval_esitmatedval_traj_sum_dssref.pdf"),
            limitsize=FALSE)
  }
}
##plotting individual ppm
for(datasourcehere in datasources){
  for(samplelevelhere in samplelevels){
    plist=vector(mode="list")
    for(compdhere in compds){
      tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")) & compd==get("compdhere"))
      if(dim(tab)[1]==0){
        next
      }
      name=paste0(compdhere,"_samplelevel",samplelevelhere,"_",datasourcehere)
      p<-ggplot(data=tab,aes(trueppm,ppm,color=matchppm))+
            geom_point(alpha=0.75)+
            xlab("true ppm")+
            ylab("traced ppm")+
            ggtitle(name)+
            guides(color=guide_legend(title="ridge_ppm"))
            theme_bw()+
            theme(plot.title=element_text(hjust=0.5))
      if(compdhere=="Acetate"){
        p=p+scale_x_continuous(limits=c(1.9,2.05))+scale_y_continuous(limits=c(1.9,2.05))
      }
      ggsave(plot=p,file=paste0(dir.r.res,name,".trueppm_esitmatedppm_traj_dssref.pdf"))
      if(compdhere %in% ppmchtable[,"compd"]){
        p=p+theme(legend.position="none",plot.title=element_blank(),
                  axis.title.x=element_blank(),axis.title.y=element_blank(),
                  axis.text.x=element_text(angle=90,hjust=1))
        plist[[compdhere]]=p
      }
    }
    ##plot a combined figure
    p<-plot_grid(plotlist=plist,labels=c(names(plist)),
            label_size=10,ncol=ncol,label_x=0.3)
    ggsave(plot=p,filename=paste0(dir.r.res,"samplelevel",samplelevelhere,"_",datasourcehere,"trueppm_esitmatedppm_traj_sum_dssref.pdf"),
            limitsize=FALSE)
  }
}

load(paste0(dir.r.res,"eval.ridtracing.RData"))
specres=readMat(paste0(dir.data,"spectral.complex.more.mat"))
ppm=specres[["ppm"]][1,]
strdata=specres[["strdata"]]
Peakofinterest=c("alanine_1"="3.76133","Serine_2"="3.91635","Glycerol_11"="3.77877","ethanol_4"="3.62311")
datasourcehere="manual"
samplelevelhere=4
nsample=50
timevec=1:nsample
threhold=0.01
tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")))
concvecstr=comptruevallist[[2]]
compose.mat.list=vector(mode="list")
compdspec_name=attr(strdata,"dimnames")[[1]]
for(peakname in names(Peakofinterest)){
  peakppm=Peakofinterest[[peakname]]
  peakind=which(tab[,"matchppm"]==peakppm)
  realppm=tab[peakind,"ppm"]
  realtime=tab[peakind,"time"]
  realintensity=tab[peakind,"intensity"]
  # realppm=realppm[realtime]
  # realintensity=realintensity[realtime]
  tempmat=as.data.frame(matrix(NA,nrow=length(realtime),ncol=3))
  colnames(tempmat)=c("time","intensity","compd")
  list.tempmat=vector(mode="list")
  tempmatsum=tempmat
  tempmatsum[,"time"]=realtime
  tempmatsum[,"compd"]="sum"
  tempmatsum[,"intensity"]=0
  for(compd in compdspec_name){
    compdind=which(compdspec_name==compd)
    concvec=concvecstr[[1]][compdind,,][[1]][1,]
    spectrum=strdata[[compdind]][,1]
    ppmind=sapply(realppm,function(x){
      which(ppm==x)
    })
    intvec=spectrum[ppmind]*concvec[realtime]
    if(max(intvec)<max(concvec)*threhold){
      next
    }
    tempmat[,"intensity"]=intvec
    tempmat[,"time"]=realtime
    tempmat[,"compd"]=compd
    tempmatsum[,"intensity"]=tempmatsum[,"intensity"]+intvec
    list.tempmat[[compd]]=tempmat
  }
  list.tempmat[["sum"]]=tempmatsum
  tempmatestimated=tempmat
  tempmatestimated[,"intensity"]=realintensity
    tempmatestimated[,"time"]=realtime
  tempmatestimated[,"compd"]="estimated"
  list.tempmat[["estimated"]]=tempmatestimated
  fullmat.peak=Reduce('rbind',list.tempmat)
  compose.mat.list[[peakname]]=fullmat.peak
  p<-ggplot(data=fullmat.peak,aes(time,intensity,color=compd))+
            geom_line()
  ggsave(plot=p,filename=paste0(dir.r.res,"samplelevel",samplelevelhere,"_",datasourcehere,"_",peakname,"deconv_scatter_plot.pdf"),
         limitsize=FALSE)
}
#calculate the related peak intensity plot the contribution > 0.01

###example plot for L-lecucine
rangelist=list("region1"=c(0.8,1.0),"region2"=c(1.5,2.0),"region3"=c(3.5,4.0))
datasourcehere="manual"
samplelevelhere=4
compdhere="Leucine"
tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")) & compd==get("compdhere"))
name=paste0(compdhere,"_samplelevel",samplelevelhere,"_",datasourcehere)
ppmregion=rep(NA,times=dim(tab)[1])
temp=sapply(names(rangelist),function(x){
  rang=rangelist[[x]]
  ind=tab[,"ppm"]>rang[1]&tab[,"ppm"]<rang[2]
  ppmregion[ind]<<-x
})
tab$ppmregion=ppmregion
tab[,"ppmregion"]=as.factor(tab[,"ppmregion"])
# p<-ggplot(data=tab,aes(time,intensity,shape=ppmregion))+
#       geom_point(alpha=0.5)+
#       xlab("time(h)")+
#       ylab("scaled_intensity")+
#       ggtitle(name)+
#       guides(color=guide_legend(title="ridge_ppm"))
#       theme_bw()+
#       theme(plot.title=element_text(hjust=0.5))
# ggsave(plot=p,file=paste0(dir.r.res,name,".timetraj.ppmregion.pdf"))
p<-ggplot(data=tab,aes(trueval,intensity,color=ppmregion))+
      geom_point(alpha=0.5)+
      xlab("true value")+
      ylab("scaled_intensity")+
      ggtitle(name)+
      guides(color=guide_legend(title="ridge_ppm"))+
      theme_bw()+
      theme(plot.title=element_text(hjust=0.5))
ggsave(plot=p,file=paste0(dir.r.res,name,".trueval_esitmatedval_traj.ppmregion.pdf"))

###exmaple plot for L-alanine
rangelist=list("region1"=c(1,3,1.6),"region2"=c(3.6,4.0))
datasourcehere="manual"
samplelevelhere=4
compdhere="alanine"
tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")) & compd==get("compdhere"))
name=paste0(compdhere,"_samplelevel",samplelevelhere,"_",datasourcehere)
ppmregion=rep(NA,times=dim(tab)[1])
temp=sapply(names(rangelist),function(x){
  rang=rangelist[[x]]
  ind=tab[,"ppm"]>rang[1]&tab[,"ppm"]<rang[2]
  ppmregion[ind]<<-x
})
tab$ppmregion=ppmregion
tab[,"ppmregion"]=as.factor(tab[,"ppmregion"])
# p<-ggplot(data=tab,aes(time,intensity,shape=ppmregion))+
#       geom_point(alpha=0.5)+
#       xlab("time(h)")+
#       ylab("scaled_intensity")+
#       ggtitle(name)+
#       guides(color=guide_legend(title="ridge_ppm"))
#       theme_bw()+
#       theme(plot.title=element_text(hjust=0.5))
# ggsave(plot=p,file=paste0(dir.r.res,name,".timetraj.ppmregion.pdf"))
p<-ggplot(data=tab,aes(trueval,intensity,color=ppmregion))+
      geom_point(alpha=0.5)+
      xlab("true value")+
      ylab("scaled_intensity")+
      ggtitle(name)+
      guides(color=guide_legend(title="ridge_ppm"))+
      theme_bw()+
      theme(plot.title=element_text(hjust=0.5))
ggsave(plot=p,file=paste0(dir.r.res,name,".trueval_esitmatedval_traj.ppmregion.pdf"))

##visualize on plotly
load(paste0(dir.r.res,"eval.ridtracing.global.dssref.RData"))
tab=filter(riddataframe,datasource==get("datasourcehere") & samplelevel==as.character(get("samplelevelhere")))
name=paste0("global_intensity_samplelevel",samplelevelhere,"_",datasourcehere)
# color=factor(ifelse(tab[,"compd"] %in% c("ethanol"),"ethanol","others"))
color=factor(tab[,"compd"])
shape=factor(tab[,"matchppm"])
p<-ggplot(data=tab,aes(trueval,intensity,color=color,shape=shape))+#
      geom_point(alpha=0.5)+
      # scale_colour_manual(values=c("ethanol"="red","others"="grey"))+
      xlab("true intensity")+
      ylab("traced intensity")+
      xlim(-100,4400)+
      ylim(-100,4400)+
      ggtitle(name)+
      # guides(color=guide_legend(title="ridge_ppm"))+
      theme_bw()+
      theme(plot.title=element_text(hjust=0.5))+
      theme(legend.position="none",
            axis.text.x=element_text(angle=90,hjust=1))
p<-ggplotly(p)

##histogram of ppm distance
file.tab=paste0(comp,"Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/result/peakmatching_simulated/manual_peak_match_differences.txt")
tab=read.table(file.tab,sep="\t",header=TRUE)
colnames(tab)=c("peaks","alpha","beta","gamma","delta")
ltab=melt(tab,id.vars=c("peaks"))
colnames(ltab)=c("peaks","data_set","value")
ltab=ltab[which(ltab[,"value"]!= -1),]
p<-ggplot(data=ltab,aes(value))+
      geom_histogram(fill="black",alpha=0.5,position="identity")+
      facet_grid(rows=vars(data_set))+
      xlab("ppm_RMSD")+
      ylab("count")+
      ggtitle("ppm_RMSD_histogram")+
      # guides(color=guide_legend(title="ridge_ppm"))
      theme_bw()+
      theme(plot.title=element_text(hjust=0.5))+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour="black"))
ggsave(plot=p,file=paste0(dir.r.res,"ppm.histogram.pdf"))
