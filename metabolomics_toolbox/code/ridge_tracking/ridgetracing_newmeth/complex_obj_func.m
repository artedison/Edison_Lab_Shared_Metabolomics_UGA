function [objstr]=complex_obj_func(mat,ppm,times,Sample,noiserang,prop_threhs,ppmpropreg,regionsele,nearbythred,dssvec)
%% this program will evaluate the complex level of spectral
%% Different spectral will be evaluated by SNR and other evaluation functions
%% The evaluation is based on ridge tracing result and so the higher the coverage of region and ridge tracing, the more precise the estimation
%% argument:
%%% mat: the input spectral matrix size(mat)=[time ppm]
%%% ppm: the ppm of the spectral, expected to be a size(ppm)=[1 length(spectra)]
%%% times: the time of each sample, expected length(time)=size(mat,1)
%%% Sample: the structure containing ridges. The manual ridge tracing result is recommended.
%%%% it is supposed to be data of one sample from ridge tracing script, which contains multiple time points
%%% noiserang: the region to calculate noise. Default [-0.4 -0.2] in ppm
%%% prop_threhs: the ratio threhold for nearby peaks. default 10
%%% ppmpropreg: the considerred spectral region in PPM. default [-1 10]
%%% regionsele: the selected region in ridge tracing step, consierred as probled region for the spectra
%%% nearbythred: the threhold to consider two peaks are close in PPM direction to calcualte Objcomdyn. default 0.05
%%% dssvec: the intensity of dss peak
%% return: objstr is a structure, including SNR, Objcomppm, Objcomdyn, Objcomrange, peakdensity
%%% SNR: signal to noise ratio. mean signal to noise level among all samples.
%%%% For each sample, SNR=mean(intensity_of_dss_peak/sd_of_selected_noise_region)
%%% Objcomppm(Shift complexity): measure ppm shift
%%%% the mean squared standard deviation of peak ppm for different peaks
%%%%% mean(sd^2(peak_ppm_each_peak))
%%%%% the peak ppm default [-1 10] is scaled to [0 1] before calculation
%%% Objcomdyn(Scale complexity): measure how dynamic are peaks in each spectral
%%%% highdiffpeakpair/(Npeak-1)
%%%% highdiffpeakpair: 1 if nearby peaks are different in intensity by prop_threhs
%%%%% 0 otherwise. Mean is calculated for all samples
%%%% nearby peaks are defined as recorded peaks within ppm differences (0.1)
%%% Objcomrange(Dynamics complexity): measure peak changes through time
%%%% the mean squared standard deviation of peak intensity proportion for different peaks
%%%%% mean(sd^2(peak_intensity_proportion_each_peak))
%%%%% the peak intensity is proported to largest intensity before calculation
%%% peakdensity: Npeak/ppmrag
%%%% Npeak: number of peak picked totally
%%%% ppmrag: the sum of ppm range for peak picking
%%%% the peak ppm default [-1 10] (assumed range) is scaled to [0 1] before calculation
%%% Attention:
%%%% baseline need to be 0 before running this function
%% Yue Wu 01/23/2019

%%%%test%%%%
% times=time';
% Sample=Sample(samp_i);
%%%%%%%%%%%%

%%parameters
if ~exist('noiserang', 'var')
  noiserang=[-0.4 -0.2];
end
if ~exist('prop_threhs', 'var')
  prop_threhs=10;
end
if ~exist('ppmpropreg', 'var')
  ppmpropreg=[-1 10];%assumed range [-1 10], converted to [0 1]
end
if ~exist('nearbythred', 'var')
  nearbythred=0.05;
end
if ~exist('dssvec', 'var')
  dssvec=[];
end
%%initialization
data=Sample.ridges;
%%%remove empty rows:
%%% might be first row or introduced manually by user while doing ridge tracing in experimental dataset
indnull=[];
for i=2:length(data)
  if length(data(i).parameters)==0
    indnull=[indnull i];
  end
end
data(indnull)=[];
ncol=length(data)-1;
Npeak=ncol;
nrow=length(times);
intensitymat=zeros([nrow ncol]);
ppmmat=zeros([nrow ncol]);
ppmcell={};
%%%reformat the intensity data in matrix and ppm data in a cell
for i=2:(ncol+1)
  tempdata=data(i).result;
  timevec=tempdata.time;
  intvec=tempdata.intensity;
  ppmvec=tempdata.ppm;
  for timeelei=1:length(times)
    timeele=times(timeelei);
    ind=find(timevec==timeele);
    if length(ind)~=0
      intensitymat(timeelei,i-1)=intvec(ind);
      ppmmat(timeelei,i-1)=ppmvec(ind);
    end
  end
  ppmcell{i-1}=tempdata.ppm;
end
minreal=min(intensitymat(:));
% intensitymat(find(intensitymat==0))=minreal;
intensitymat=intensitymat-minreal;%%shift spectral peak matrix to remove negative value
noisereg=matchPPMs(noiserang,ppm);
noisemat=mat(:,noisereg(1):noisereg(2));
%% signal to noise ratio
% refpeakind=find(ppmmeanvec>rangdss(1)&ppmmeanvec<rangdss(2));
% intensitymat(:,refpeakind)
if length(dssvec)~=0
  SNR=mean(dssvec./std(noisemat,0,2));
else
  SNR=max(intensitymat(:))/std(noisemat(:));
end
%% measure ppm shift
ppmvarvec=[];
for i=1:ncol
  scalppm=(ppmcell{i}-ppmpropreg(1))/(ppmpropreg(2)-ppmpropreg(1));
  scalppm_sq=sum((scalppm-mean(scalppm)).^2)/length(scalppm);
  ppmvarvec=[ppmvarvec scalppm_sq];
end
Objcomppm=mean(ppmvarvec);
%% measure how dynamic are peaks in each spectral
ppmmeanvec=mean(ppmmat,1);
[ppmmeanvecsort ppmsortind]=sort(ppmmeanvec);
intensitymatsort=intensitymat(:,ppmsortind);
dynvec=[];
inddiff=ppmmeanvecsort(2:ncol)-ppmmeanvecsort(1:(ncol-1))<nearbythred;
npeakpair=length(find(inddiff));
for i=1:nrow
  specpeakvec=intensitymatsort(i,:);
  indhighdiffpeakpair=abs(log10(specpeakvec(2:ncol)./specpeakvec(1:(ncol-1))))>log10(prop_threhs);
  dynvec=[dynvec length(find(inddiff&indhighdiffpeakpair))];
end
Objcomdyn=mean(dynvec)/(npeakpair);
%% measure peak changes through time
peakprop_sdvec=[];
for i=1:ncol
  intensivec=intensitymat(:,i);
  propintpeaks=intensivec/max(intensivec);
  propintpeaks_sq=sum((propintpeaks-mean(propintpeaks)).^2)/length(propintpeaks);
  peakprop_sdvec=[peakprop_sdvec propintpeaks_sq];
end
Objcomrange=mean(peakprop_sdvec);
%% peak intensity
%%% use sum to ppm region to approximate total region without considering overlpped region
ppmrag=0;
nrowregion=size(regionsele,1);
for i=1:nrowregion
  rag=regionsele(i,:);
  ppmrag=ppmrag+rag(2)-rag(1);
end
% peakdensity=Npeak/ppmrag*(ppmpropreg(2)-ppmpropreg(1));
peakdensity=Npeak/ppmrag;
%%return
objstr.SNR=SNR;
objstr.Objcomppm=Objcomppm;
objstr.Objcomdyn=Objcomdyn;
objstr.Objcomrange=Objcomrange;
objstr.peakdensity=peakdensity;
