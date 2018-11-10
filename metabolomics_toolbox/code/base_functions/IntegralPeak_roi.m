function Auc=IntegralPeak_roi(X,ppm,ppmLeft,ppmRight)

%% calculate the area under curve for regions of interest in spectra
% Input
% X: stack of 1D spectra
% ppm: chemcal shift vector 
% ppmLeft=[2.31, 4.32];
% ppmRight=[2.50, 4.54];
% Output
% Auc: the area under curve for ROIs
%%

%if exist('Y')==0
%    Y=zeros(size(X,1),1);
%    Ycolor=ones(size(X,1),1);
%else
%Ycolor=ceil(([(Y-mean(Y))/(2.01*max(abs(Y-mean(Y))))]+.5)*100);
%end
%cmap=jet(100);

%Auc=zeros(size(X));
for k=1:length(ppmLeft)
    [~,idxL]=min(abs(ppmLeft(k)-ppm));
    [~,idxR]=min(abs(ppmRight(k)-ppm));

    if idxL>idxR
        Xseg=X(:,idxR:idxL);
    else
        Xseg=X(:,idxL:idxR);
    end
    Auc(:,k)=trapz(Xseg,2);
end

%%
%figure;
%hold on; 
%for k=1:size(X,1)
%    stem(ppm,Auc(k,:),'fill','Color',cmap(Ycolor(k),:));
%end
%set(gca,'XDir','reverse')

