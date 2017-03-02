function out = overlayImagesc(x1,x2,y1,y2,z1,z2)

if 0
    % Example
    x1=linspace(1.2, 100.4, 100);
    x2=linspace(10.1, 140.7, 100);
    y1=linspace(-11.1, 100.7, 100);
    y2=linspace(-20.6, 90.3, 100);
    z1=zeros(100);
    z2=zeros(100);
    %Put a box on value 50 (not index 50)
    z1(53:57,47:52)=1;
    z2(62:66,28:33)=1;
    figure,h1=subplot(1,2,1);, imagesc(x1,y1,z1), h2=subplot(1,2,2), imagesc(x2,y2,z2)
    linkaxes([h1,h2])
end
%%
xMin=min([x1,x2]);
xMax=max([x1,x2]);
yMin=min([y1,y2]);
yMax=max([y1,y2]);
numPointsX=min(length(x1),length(x2));
numPointsY=min(length(y1),length(y2));
numPoints = min([numPointsX,numPointsY]);
x=linspace(xMin,xMax,numPoints);
y=linspace(yMin,yMax,numPoints);

z1n=zeros(size(z1,1),length(y));
for i=1:size(z1,1)
    z1n(i,:) = interp1(y1, z1(i,:), y);
end
z1=z1n;
clear z1n
for i=1:size(z1,2)
    z1n(:,i) = interp1(x1, z1(:,i), x);
end

z2n=zeros(size(z2,1),length(y));
for i=1:size(z2,1)
    z2n(i,:) = interp1(y2, z2(i,:), y);
end
z2=z2n;
z2n=zeros(length(x),size(z2,2));
for i=1:size(z2,2)
    z2n(:,i) = interp1(x2, z2(:,i), x);
end

%figure,imagesc(x,y,z2n'+z1n)

out.ppmH=x;
out.ppmX=y;
out.z1=z2n';
out.z2=z1n;