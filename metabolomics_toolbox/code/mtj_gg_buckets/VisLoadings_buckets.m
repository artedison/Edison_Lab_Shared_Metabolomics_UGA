function VisLoadings_buckets(PCA,component,X,ppm,bucketBounds)

%% VisLoadings_buckets

% Author: MTJ
% Version: 0.1
% Date: 2020
%
% Description:
%
%   Plots buckets with shading representing colorscale that corresponds to
%   the loadings from the provided model (PCA or PLS). 
%
% Inputs:
%
%     Required Arguments:
%         PCA           Model from nipalsPCA or one of the PLS functions.
%                       Must include a loadings field.
%         component     component number to use for loadings plot
%         X             spectral matrix
%         ppm           ppm vector
%         buckets       bucket bounds 
%
%     Optional Name,Value pair Arguments: 
% 
%         'addBox'      any value works; if set to anything, then a box is drawn around the legend
%         'numColumns'  must be a whole number. Number of legend columns you want
%         'textSize'    number defining text size in legend
%         'lineWidths'  width of the lines (to match your spectra)
%
% Output:
%       
%       Spectral plot with bucket bounds drawn and bucket colors
%       corresponding to the loadings values. 
%
% Usage: 
%
%         VisLoadings_buckets(model,1,matrix,ppm,bucketBounds)
%                 
% MTJ and GG 2020

    % Several components as bars below would be nice
        %colors = [0.5686    0.0392    0.0392;1.0000    1.0000    1.0000;0  0  0.5686];
        
    % Generate the colormap
            
        v = PCA.loadings(component,:);
        rng = max([min(abs(v)),max(abs(v))]);
        [cinfo] = customColormap(v,'setPoints',[-rng,0,rng]);
        
    % Plot the bins on the data
        
        figure, hold on
            highlightROIs(bucketBounds',max(X(:)),'colorVect',cinfo.rgb,'edgeColor','none','transparency',1)
            plot(ppm,X)
            set(gca,'XDir','reverse')
            title({['Loadings for Component',num2str(component)];...
                   ['(',num2str(round(PCA.variance(component)*100,1)),'% of total variance)']})
            xlabel('Chemical Shift (ppm)')
            ylabel('Signal Intensity')
            yticks(gca,[])

            colormap(cinfo.cmap);
            t=colorbar;
            set(t,'Ticks',(cinfo.setPts-min(cinfo.setPts))/range(cinfo.setPts))
            set(t,'TickLabels',num2cell(round(cinfo.setPts,2)))
            set(get(t,'ylabel'),'String', 'Loadings')
            
            set(gca,'FontSize',20)
%             set(gca,'XDir','rev');
%             hold off
    
    %% Plot the lines with the right colors
    
        % Plot a multicolored line
            % Construct the colors list
%                 
%                 colors = zeros(length(ppm),3); % preallocate length(ppm),3
%                 y = fillRegions(buckets);
%             
%             z = zeros(size(x));
%                 col = rgb;  % This is the color, vary with x in this case.
%                 surface([ppm;ppm],[y;y],[z;z],[col;col],...
%                         'facecol','no',...
%                         'edgecol','interp',...
%                         'linew',lineWidth);

end