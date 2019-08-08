function [vec ppmnewshift]=ppmshift(spectral,baseph,acidph,pka,pH,ppm,pHmeasu)
%% this function will calculate new ppm from pH
%% it will use the equation:
%%% ppm=(ppmacidph+ppmbaseph*10^(pH-pka))/(1+10^(pH-pka));
%% to calculate ppm shift in this pH compared with those in gissmo measured pH to shift the spectra.
%%% the shift is done by move the part moving out of one bounds to another bounds
%% ref Modelling the acid/base 1H NMR chemical shift limits of metabolites in human urine
%% assumption: all peak shift are within ppm range/the baseline is random noise
%% argument:
%%% spectral: the reference spetra
%%% baseph: ppm of the complete base form
%%% acidph: ppm of the complete acid form
%%% pka: pka of the compound, only two forms of compound (L and HL) is assumed.
%%% pH: the current pH
%%% ppm: the ppm vector
%%% pHmeasu: the pH where gissmo spectra is measured default 7.4
%% return: vec, the shifted spectral
%%  ppmnewshift: the ppm shift for the peak under corresponding pH
%% by Yue Wu 11/25/2018

%%test%%%%%%%%%%%%
% baseph=parahere(3);
% acidph=parahere(2);
% pka=parahere(1);
% pH=seqpH(i);
%%%%%%%%%%%%%%%%%%
if ~exist('pHmeasu', 'var')
  pHmeasu=7.4;
end
ppm1=(acidph+baseph*10^(pHmeasu-pka))/(1+10^(pHmeasu-pka));
ppm2=(acidph+baseph*10^(pH-pka))/(1+10^(pH-pka));
ppmnewshift=ppm2-ppm1;
ind=matchPPMs([ppm1 ppm2],ppm);
deltppm=ind(2)-ind(1);
len=length(spectral);
vec=[];
if deltppm>0
  vec([(deltppm+1):len 1:deltppm])=spectral;
  vec=vec';%% to unify the dimension of vectors
else
  deltppm= -deltppm;
  vec=spectral([(deltppm+1):len 1:deltppm]);
end
