function XREF=Correct_Ref(X,ppm)
% Correct the reference value for the spectra
% Arguments:
% X     Data matrix of spectra
% ppm   Chemical shift vector corresponding to X
% Return Values:
% XREF  Spectral matrix with corrected references.

 [~,k0]=min(abs(ppm));
 [~,k1]=min(abs(ppm+0.5));
for i=1:size(X,1)
 XDSS=X(i,k1:k0);
 [maxtab,~]=peakdet(XDSS,1e5);
 [~,idx]=max(maxtab(:,2));
 kdss=k0-k1-maxtab(idx,1);
 XREF(i,:)=circshift(X(i,:),[0,kdss]);
end
