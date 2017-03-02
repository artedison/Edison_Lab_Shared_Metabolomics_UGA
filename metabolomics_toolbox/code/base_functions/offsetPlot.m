
function offsetPlot(ppm, X, multPer, color)
if nargin==2
    multPer=.01;
end

maxdiff=max(sum(X));
mult=multPer*maxdiff;

a=repmat([0:mult:mult*size(X,1)-1]',1,size(X,2));
X2=X+a;
colors=blue(size(X,1)+1);
colors=colors(1:end-1,:);

if nargin==4
    figure, hold
    for i=1:size(X,1)
        plot(ppm,X2(i,:),'Color',colors(i,:))
    end
else
    plot(ppm,X2)
end
set(gca,'xdir','rev')