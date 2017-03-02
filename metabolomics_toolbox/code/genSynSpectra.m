function [ppmX,X2]=genSynSpectra(nC,v,nR,noise)
if nargin<4
    noise=0.1;
end
if nargin<3
    nR=10;
end
if nargin<2
    v=0.1;
end
if nargin<1
    nC=5;
end

%List all ft files
cd('E:\Chaevien_Shared\C13\Standards_C\proc_file')
cd('E:\Chaevien_Shared\C13\Standards\proc_file\1H')
d=dir;
d={d.name};
d=d(strfindin('.ft',d));
ftlist = cellfun(@(x) ([cd filesep x]),d,'uniformoutput',0);
%Pick nC random ft files.
idx=randperm(length(ftlist));
ftlist=ftlist(idx(1:nC));
clear spectra
for i=1:length(ftlist)
    spectra(i)=pipe2matlab(ftlist{i});
end


%% Reference
spectra = ref_spectra_manual(spectra);
[X,ppmX,titles]=Setup1D(spectra);

%% remove region
X=remove_region(X,ppmX,4.7,5.05);
[X,ppmX]=remove_ends(X,ppmX,-.05,10);
%% Bring spectra to median of 0 (shift up or down)
X=bsxfun(@minus,X,median(X,2));

%%
i=5;
[A,S]=showBaseline(X(i,:))
showBaseline(X(i,:),ppmX,A,S)
%% normalize
XN=normalize(X,ppmX,'integral',[-.02,0.02]);

%% For each row, generate nR more with variance v. then add em all together
% The variance is applied to each compound. So each peak for a given
% compound will have perfect correlation/covariance
X2=zeros(nR,size(X,2));
for i=1:nC
    c=bsxfun(@plus,(X(i,:)'*normrnd(0,v,[1,nR]))',X(i,:));
    X2=X2+c;
end

%% Add random noise to every index at every sample. Peaks in a single compound now
% won't be perfect. Note: This will also mess up peak shapes because its
% being applied to every point, not every peak.
X2=X2+(X2.*normrnd(0,noise,size(X2)));



