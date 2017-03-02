function XALg=HATS(X,XNoise,ppm1,ppm2,label,maxshift)

% XALg=HATS(X,XNoise,ppm1,ppm2,label,maxshift)
%
% Alignment using HATS method. Robinette et al. Analytical Chemistry,
% 2011. 83(5):1649-1657.  Applies segment2D if no label matrix is provided
%
% Arguments:
% 
% X            Data matrix of spectra
% XNoise       Local Noise matrix
% ppm1         Chemical shift vector of F2
% ppm2         Chemical shift vector of F1
% label        label matrix from segment2D.m (optional)
% maxshift     Maximum shift if F2 and F1 dimensions [F2,F1] (optional)
%
% Return Values:
% XALg         Aligned spectral matrix
%
% Steven Robinette


%% Detect Image Processing Toolbox

if isempty(ver('images'))==0
    disp('You do not have the Image Processing Toolbox.  HATS will run slightly slower.')
end

XALg=X;
XALg(isnan(XALg))=0;

%% Package spectral set into cell array "align"

for ind1=1:size(X,3)
    align{ind1}=ind1;
end

%% Calculate correlation matrix of spectral bins

if exist('maxshift')~=1
    if abs(max(ppm1)-max(ppm2))<10
        maxshift=[0.1 0.1];
    else
        maxshift=[1 0.1];
    end 
end

thresh=1000;
if exist('label')~=1
    label=segment2D(X,XNoise,ppm1,ppm2,thresh);
end

binmat=zeros(size(X,3),max(max(label)));
label3d=zeros(size(X));
for k=1:size(X,3);
label3d(:,:,k)=label;
end

for k=1:max(max(label))
indices=find(label3d==k);
intensitymatrix=reshape(X(indices),length(indices)/size(X,3),size(X,3));
binmat(:,k)=sum(abs(intensitymatrix));
end


%% Create spectral guide tree

m=linkage(binmat,'weighted','spearman');

figure, h=dendrogram(m,0,'ORIENTATION','left');
title('HATS Guide Tree')

%% Align up dendrogram
tic

for aligndex=1:size(m,1)
    
    num1=m(aligndex,1);
    num2=m(aligndex,2);
   
    XALg(:,:,[align{num1},align{num2}])=star_align2D(XALg(:,:,[align{num1},align{num2}]),label,ppm1,ppm2,maxshift);
    sum(sum(sum([XALg(:,:,align{num1})-X(:,:,align{num1})])))
    
    align{ind1+1}=[align{num1},align{num2}];
    
    ind1=ind1+1;
    disp(['Completed ',num2str(aligndex),' of ',num2str(size(m,1))]);
    toc
end
toc

