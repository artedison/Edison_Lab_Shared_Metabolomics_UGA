
function spectra=ref_chae(spectra,thresh)
%% Pick all peaks
if ~exist('thresh','var')
    thresh = 0.01;
end
for i=1:length(spectra)
    a=0;
    peaks=peakpick(spectra(i).real,30,thresh*max(spectra(i).real));
    
    subplot(length(spectra),1,i); plot(spectra(i).ppm,spectra(i).real,'k')
    hold
    plot(spectra(i).ppm(peaks),spectra(i).real(peaks),'ro')
    title(spectra(i).Title);
    set(gca,'XDir','rev')
    
    ax(i)=subplot(length(spectra),1,i);
    % find the peak closest to 0 and set it equal to zero
    
    [~, index] = min(abs(spectra(i).ppm(peaks) - (-0.65)));
    
    TSP_int(i)=spectra(i).real(peaks(index));
    
    TSP_ppm(i)=spectra(i).ppm(peaks(index));
    
    if TSP_ppm(i)<0;
        spectra(i).ppm=spectra(i).ppm+abs(TSP_ppm(i));
    else
        spectra(i).ppm=spectra(i).ppm-abs(TSP_ppm(i));
    end
    
    clear peaks index
    
end
linkaxes(ax(1,:),'x')

%% Isolate the TSP peak in preparation for integration
% for i=1:length(spectra)
%     
%     b=0.02;
%     [~,index(i,1)]=min(abs(spectra(i).ppm - b));
%     TSP_ppm(i,1)=spectra(i).ppm(index(i,1));
%     
%     [~,index(i,2)]=min(abs(spectra(i).ppm + b));
%     TSP_ppm(i,2)=spectra(i).ppm(index(i,2));
%     
% end
%% Plot

% for i=1:length(spectra)
%     Reference={'Reference Peak'};
%     figure, plot(spectra(i).ppm(index(i,1):index(i,2)),spectra(i).real(index(i,1):index(i,2)));
%     title([spectra(i).Title,Reference]);
%     
% end

%% Isolate the TSP from the spectra then zero the TSP in the Spectra
% 
% for i=1:length(spectra)
%     ppm(:,i)=spectra(i).ppm(index(i,1):index(i,2));
%     real(i,:)=spectra(i).real(index(i,1):index(i,2));
% end
% 
% % zero in spectra
% for i=1:length(spectra)
%     spectra(i).real(index(i,1):index(i,2))=0;
% end
% 

%% Divide the TSP by its area to give an area of 1 for the TSP
% 
% for i=1:length(spectra)
%     
%     % Uses peakfit to find the area under the TSP peak
%     figure, [FitResults]=peakfit([ppm(:,i) real(i,:)'],0,0.1,1,2,0,5);
%     
%     TSP_area(i,1)=FitResults(1,5);
%     
%     % Divide the TPS peak area by the itself to give a peak area of approximately 1, TSP_area(:,1)=old peak area and TSP_area(:,2)=new peak area
%     real(i,:)=real(i,:)./TSP_area(i,1);
%     spectra(i).real=spectra(i).real./TSP_area(i,1);
%     
%     figure, [FitResults]=peakfit([ppm(:,i) real(i,:)'],0,0.1,1,2,0,5);
%     
%     TSP_area(i,2)=FitResults(1,5);
% end
% 
% % Put TSP back into the spectra
% for i=1:length(spectra)
%     spectra(i).real(index(i,1):index(i,2))=real(i,:);
% end
 

end


%%

 