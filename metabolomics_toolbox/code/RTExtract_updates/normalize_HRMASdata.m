function [matrixN,ROInorm,factors] = normalize_HRMASdata(matrix,ppm,method,ROInorm)
%% Normalize based on a peak

%% Select the peak to normalize on


if isempty(ROInorm)
    h = figure, plot(ppm,matrix)
        set(gca,'XDir','reverse')
        title('Select Peak to Normalize On')
        xlabel('Chemical Shift (ppm)')
        ylabel('Signal Intensity')
        ROInorm = [];
        whichLine()

    selectROIsFromFigure();
    D=get(gcf,'Children'); %get the handle of the line object
        patches = findall(D,'Type','Patch');  % make a list of all patch boxes in the figure
        heights = zeros(1,length(patches));
        ROInorm = zeros(2,length(patches));
            for i = 1:length(patches)       % for each patch box, get the vertices and store
                verts = patches(i).Vertices;
                ROInorm(:,i) = [verts(1);verts(3)];
                heights(i) = abs(verts(3)-verts(5));
            end
        close(h)
   
end       
%% Do the normalization depending on the method
    newROIinds = matchPPMs(ROInorm,ppm);
    switch method
        case 'maximum'
            factors = max(matrix(:,newROIinds(1):newROIinds(2)),[],2);

        case 'integral'
            factors = sum(matrix(:,newROIinds(1):newROIinds(2)),2);
    end

    matrixN = matrix./repmat(factors,1,size(matrix,2));
    
end