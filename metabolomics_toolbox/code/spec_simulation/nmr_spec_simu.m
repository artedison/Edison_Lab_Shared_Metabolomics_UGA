function [spec_mat] = nmr_spec_simu(spec_args)
% the function simulate spectra based concentration matrix (sample X compound) and gissmo spectra
% both white noise and linebroadening will be introduced
% Arguments: spec_args:
%            spec_args.ppm: the ppm vector
%            spec_args.sampleindvec: sample index vector. EX. [1 1 1 2 2 2] two groups and each with 3 samples
%            spec_args.compd_vec: the simulated compound list. Must be provided.
%            spec_args.tabinfor: the information table for downoaded libarary. Must beprovided
%            spec_args.conc_mat: the concentration matrix from function spec_conc_simu()
%            spec_args.strdata: the gissmo reference spectra.
%            spec_args.sigma_noise: white noise on the whole spectra. default 15.
%            spec_args.lambda: control line broadening. default -0.00035. the exponential window function for fid.
%            spec_args.totalnoiseregion: the region of no signal (to shift spectra). default [-0.5 -0.2]
% Return: spec_mat: the simulated spectra for multiple samples
% YUE WU 12112019

if ~isfield(spec_args,'ppm')
  error('please provide ppm vector');
end
if ~isfield(spec_args,'sampleindvec')
  error('please provide sample vector');
end
if ~isfield(spec_args,'compd_vec')
  error('please provide compound list');
end
if ~isfield(spec_args,'tabinfor')
  error('please provide library information');
end
if ~isfield(spec_args,'conc_mat')
  error('please provide concentration matrix');
end
if ~isfield(spec_args,'strdata')
  error('please provide reference spectra');
end
if ~isfield(spec_args,'sigma_noise')
  spec_args.sigma_noise=15;
end
if ~isfield(spec_args,'lambda')
  spec_args.lambda=-0.00035;
end
if ~isfield(spec_args,'totalnoiseregion')
  spec_args.totalnoiseregion=[-0.5 -0.2];
end

ppm=spec_args.ppm;
sampleindvec=spec_args.sampleindvec;
compd_vec=spec_args.compd_vec;
tabinfor=spec_args.tabinfor;
conc_mat=spec_args.conc_mat;
strdata=spec_args.strdata;
sigma_noise=spec_args.sigma_noise;
lambda=spec_args.lambda;
totalnoiseregion=spec_args.totalnoiseregion;

compdvec_all=table2cell(tabinfor(:,1));
noiserang=matchPPMs(totalnoiseregion,ppm);
noiseppmind=noiserang(1):noiserang(2);
spec_mat=zeros(length(sampleindvec),length(ppm));
for samplei=1:length(sampleindvec)
  matvechere=zeros(1,length(ppm));
  for compdi=1:length(compd_vec)
    compd=compd_vec{compdi};
    compdind=find(strcmp(compdvec_all,compd));
    inforvec=table2cell(tabinfor(compdind(1),:));
    ratio=conc_mat(samplei,compdi);
    name=[inforvec{2} '_' inforvec{3}];
    compd_spec=strdata.(name);
    matvechere=matvechere+compd_spec'*ratio;
  end
  matvechere=matvechere+random('Normal',0,sigma_noise,size(matvechere,1),size(matvechere,2));
  Y=hilbert(matvechere);
  X1=ifft(Y);
  timesresolu=1:1:length(matvechere);% the first point should be fixed to 1. normalization
  X1=X1.*exp(lambda*(length(matvechere)-timesresolu));
  y1=fft(X1);
  new_matvechere=real(y1);
  %%shift the intercept of spectral
  new_matvechere=new_matvechere-mean(new_matvechere(noiseppmind));
  spec_mat(samplei,:)=new_matvechere;
end
