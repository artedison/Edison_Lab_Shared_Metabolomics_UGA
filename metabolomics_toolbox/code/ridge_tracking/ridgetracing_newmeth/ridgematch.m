function [res]=ridgematch(ppm,ridarray)
%% this function is for match ridges to peaks list of multiple compounds
%% for each peak(ppm) of the compound, the nearest ridges will be chosen.
%% relative intensity between differnet peaks for the same compound is not used as information here yet
%% the mse is normalized by length of ridge
%% argument:
%%% peak: the ppm to match
%%% ridarray: the ridges structure produced from the function 'ridgetrace_power2_ext'
%% return: res a structure that contains the candidate ridge and corresponding mse
%%% mse, the minimum squared error between ridges and the expected peak location
%%%%mse: mean squared error
%%%%sse: sum squared error
%% YW 12/1/2018

% %%%%%test%%%%%%%
% ppm=comppeaklist(t,1);
% ppm=comppeaklist;
% ridarray=rids;
% ppm=ppmcandidate;ridarray=rids;
% %%%%%%%%%%%%%%%%

sqdist=[];
realmse=[];
ssevec=[];
for j=2:length(ridarray)%%first row of the ridge array is empty as default
  ppmvec=ridarray(j).result.ppm;
  ppmvec=ppmvec(:);
  %%the result from ridge tracing might be any order, though correspondingly.
  %%this step will make sure the two ppm vector have same order
  if length(ppm)~=1
    indreform=ridarray(j).result.rowind;
    ppmreform=ppm(indreform);
    ppmreform=ppmreform(:);
  else
    ppmreform=ppm(:);
  end
  len=length(ppmvec);
  sse=sum((ppmvec-ppmreform).^2);%sum squared error
  sqdist=[sqdist sse/len];%%scaled sse
  ssevec=[ssevec sse];%%sse
end
[scalmse ind]=min(sqdist);
res=struct();
res.candidate=ridarray(ind+1);
res.scalmse=scalmse;
res.sse=ssevec(ind);
res.ind=ind;
end
