function [corr,covar]=chaeSTOCSY(target,X,ppm)

% STOCSY(target,X,ppm)
% 
% Plots correlation/covariance projection of NMR spectrum to target
% chemical shift or response vector

% Arguments:
% 
% target       Chemical shift of driver peak or response vector with length equal to size(X,1)
% X            Data matrix of spectra
% ppm          Chemical shift vector corresponding to X
%

if length(target)==1
    [h,k]=min(abs(ppm-target));
    target_vect=X(:,k);
else
    target_vect=target;
end
    
if size(target_vect,1)~=size(X,1) && size(target_vect,2)==size(X,1)
    target_vect=target_vect';
end

corr=(zscore(target_vect')*zscore(X))./(size(X,1)-1);   
covar=(target_vect-mean(target_vect))'*(X-repmat(mean(X),size(X,1),1))./(size(X,1)-1);
% 
% cmap=jet(100);

% lines=NaN(size(cmap,1),size(corr,2));
% ind=1;
% for k=-1:2/size(cmap,1):.99
% lines(ind,find(corr>k))=covar(find(corr>k));
% ind=ind+1;
% end
% 
% figure, plot(ppm,lines(1,:),'Color',cmap(1,:))
% hold on
% for k=2:size(cmap,1)
%     plot(ppm,lines(k,:),'Color',cmap(k,:));
% end
% set(gca,'XDir','rev')
% xlabel('Chemical Shift (ppm)')
% if length(target==1)
%     ylabel(['Covariance with signal at ',num2str(target),' ppm'])
% else
%     ylabel('Covariance with Y vector')
% end
% 
% t=colorbar;
% if length(target==1)
%     set(get(t,'ylabel'),'String', ['Correlation with signal at ',num2str(target),' ppm']);
% else
%     set(get(t,'ylabel'),'String', ['Correlation with Y vector']);
% end
% 
% caxis([-1 1])