function overlayNMR2D(spectra,thresh)
% accepts a structure array containing real, ppm1, ppm2, Title

if exist('thresh')==0
    thresh=5;
end

if length(spectra)>length(thresh)
    thresh=thresh(1);
    thresh=repmat(thresh,length(spectra),1);
end
    
    
for k=1:length(spectra)
    for y=1
     
        range=3;
        levels=10;
 
        subplot(2,length(spectra),k)  
        colors=('rbgkmyrbgkmy');
            
        vector=(2.^[-1*range:(range-(-1*range))/(levels-1):range])*(thresh(k)*std(std(spectra(y,k).real)));
        h=contour(spectra(y,k).ppm1,spectra(y,k).ppm2,spectra(y,k).real,vector,'EdgeColor',colors(k));
        set(gca,'XDir','rev')
        set(gca,'YDir','rev')
        xlabel('F2 (ppm)');
        ylabel('F1 (ppm)');
        title(spectra(y,k).Title(1,:))
        
        ax(k)=subplot(2,length(spectra),k);
    end
end

for k=1:length(spectra)
    for y=1
     
        range=3;
        levels=10;
        subplot(2,length(spectra),[length(spectra)+1:2*length(spectra)])  
        colors=('rbgkmyrbgkmy');
            
        vector=(2.^[-1*range:(range-(-1*range))/(levels-1):range])*(thresh(k)*std(std(spectra(y,k).real)));
        h=contour(spectra(y,k).ppm1,spectra(y,k).ppm2,spectra(y,k).real,vector,'EdgeColor',colors(k));
        set(gca,'XDir','rev')
        set(gca,'YDir','rev')
        xlabel('F2 (ppm)');
        ylabel('F1 (ppm)');
        title('Overlayed Spectra')
        ax(length(ax)+1)= subplot(2,length(spectra),[length(spectra)+1:2*length(spectra)]);  
    end
    hold on
    ishold;
    
end

linkaxes(ax(1,:),'xy')