function [flag_slope,flag_residuals]=autoflag_reg(run_order,X,XTitles,metaboliteList)
%%

% X is a data matrix of 2 dimenstions and Y is a vector
%%

Y=run_order
flag_residuals=XTitles


if isstr(metaboliteList)==0
    metaboliteList1=num2str(metaboliteList)
    metaboliteList2=strsplit(metaboliteList1)
else
end

flag_slope=metaboliteList2'


for i=1:length(X(1,:))
    fit=fitlm(Y,X(:,i))
    pvalue=fit.Coefficients{2,4}
     for j=pvalue<0.05
       flag_slope(i,2)={j}
     end
    R=fit.Residuals{:,4}
    for k=abs(R(:,1))>1
        flag_residuals=horzcat(flag_residuals,k)
    end
h=figure 
hold on
subplot(1,2,1) 
hold on
scatter(Y((flag_residuals(:,i+1)==0)),X(((flag_residuals(:,i+1)==0)),i),'b')
scatter(Y((flag_residuals(:,i+1)==1)),X(((flag_residuals(:,i+1)==1)),i),'r')
hold off
xlabel('Run Order')
ylabel('AU')
title('Plotted against Run Order - Flagged Residuals in Red','FontSize', 10)
subplot(1,2,2)
plotResiduals(fit,'probability','ResidualType','Standardized')
hold off

name=metaboliteList2{i}
print (h, '-dpdf',name)
print (h, '-depsc',name)

end




