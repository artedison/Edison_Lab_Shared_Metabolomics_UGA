function binmat=bin2D(X,label)

% binmat=bin2D(X,XNoise)
%
% Automatic binning of 3D spectral matrix X using spectral segments in
% label matrix
%
% Arguments:
% X                    Spectral matrix
% label                Spectral segmentation matrix
%
% Return Values:
% binmat               Matrix of bin intensities
%
% Steven Robinette

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

