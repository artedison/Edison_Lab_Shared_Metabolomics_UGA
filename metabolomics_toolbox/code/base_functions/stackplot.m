function stackplot(X,XNoise,ppm1,ppm2,Y,thresh)

% stackplot(X,XNoise,ppm1,ppm2,Y,thresh)
% 
% Plots stacked contour plot of X colored by Y (if known).
%
% Arguments:
% 
% X            Data matrix of spectra
% XNoise	   Calculated noise matrix
% ppm1         Chemical shift vector of F2
% ppm2         Chemical shift vector of F1
% Y            if response vector Y is known, use it to color the contour
%              plots (optional)
% thresh       Intensity threshold- can be scalar or vector; if vector
%              same number of elements as size(X,3)

if exist('thresh')~=1
    for k=1:size(X,3)
        test=X(:,:,k)./XNoise(:,:,k);
        threshvect(k)=mean(test(find(test~=0)))+5*std(test(find(test~=0)));
    end
elseif length(thresh)>1
    threshvect=thresh;
    if size(threshvect,1)>size(threshvect,2)
        threshvect=threshvect';
    end
    clear thresh
else
    threshvect=repmat(thresh,1,size(X,3));
end

if exist('Y')==0
    [junk,Ycolor]=sort(rand(size(X,3),1));
else
    Ycolor=ceil(([(Y-mean(Y))/(2.01*max(abs(Y-mean(Y))))]+.5)*100);
end
cmap=jet(100);


range=3;
levels=10;
vector=(2.^[-1*range:(range-(-1*range))/(levels-1):range])*threshvect(1);

figure
h=contour(ppm1,ppm2,X(:,:,1)./XNoise(:,:,1),vector,'EdgeColor',cmap(Ycolor(1),:));
set(gca,'XDir','rev')
set(gca,'YDir','rev')
for k=2:size(X,3)
    thresh=threshvect(k);
    vector=(2.^[-1*range:(range-(-1*range))/(levels-1):range])*thresh;
    hold on
contour(ppm1,ppm2,X(:,:,k)./XNoise(:,:,k),vector,'EdgeColor',cmap(Ycolor(k),:))
disp(k)
end