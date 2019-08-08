function [matreturn ppmrealvec]=poly_mov_spec_syn(refspectra,ppm,refpeakppm,sampsize,timevec,ppmvec,npower,direction,conchigh)
%% this function will move a the peak of a single peak reference spectral according to a npower polynomial function
%% the polynomial will be ppmvec=f(timevec) with power npower
%% argument:
%%% refspectra: the reference spectral to base on
%%% ppm: the ppm vector of the reference spectra
%%% refpeakppm: the peak position of the reference spectra
%%% sampsize: the sample size to synthize
%%% timevec: the time vector for polynomial
%%% ppmvec: the ppm vector for polynomial
%%% npower: the power of the polynomial
%%% direction: the direction of changes of peak intensity. 1 increase, 0 no change, -1 decrease
%%% conchigh: the highest intensity of the moving peak. It will be the intenstiy of peaks if drection=0
%% return: a spectral matrix with each row different time points(samples) and each colume different ppm
%%%% yue wu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%test%%%%%%
%%%%%%%%%%%%%%%

%% this part synthize a spectral named as unknown
p=polyfit(timevec,ppmvec,npower);
ppmrealvec=polyval(p,1:sampsize);
matreturn=zeros([sampsize length(ppm)]);
%%different chemical shift of spectra
for i=1:sampsize
  ppm2=ppmrealvec(i);
  ind=matchPPMs([refpeakppm ppm2],ppm);
  deltppm=ind(2)-ind(1);
  len=length(refspectra);
  if deltppm>0
    vec([(deltppm+1):len 1:deltppm])=refspectra;
    vec=vec';
  else
    deltppm= -deltppm;
    vec=refspectra([(deltppm+1):len 1:deltppm]);
  end
  matreturn(i,:)=vec;
end
%%different intensity of the spectra
if direction==1
  ratio=linspace(0,conchigh,sampsize);
elseif direction== -1
  ratio=linspace(conchigh,0,sampsize);
else%%direction==0
  ratio=repmat(conchigh,[1 sampsize]);
end
matreturn=matreturn.*ratio';
