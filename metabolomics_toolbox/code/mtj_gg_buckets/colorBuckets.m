function [ax] = colorBuckets(X,ppm,buckets,bucketColors,titleStr,colorbarStr,transparency,lineColors,lineLabels)
%% colorBuckets

% Author: MTJ
% Version: 0.1
% Date: 2020
%
% Description
%
%   Plots buckets using colors provided in the cinfo structure (e.g. from  
%   customColormap()). Useful for plotting loadings or ttest statistics or 
%   correlations on bucketed data.
%
% Inputs:
%
%     Required Arguments:
%         'X'                   spectral matrix
%         'ppm'                 ppm vector
%         'buckets'             bucket boundaries (n x 2)  
%         'cinfo'               cinfo data structure, e.g. from customColormap()
%         'titleStr'            plot title string (use '' for no title)
%         'colorbarStr'         label for colorbar (e.g. 'Correlation')
%         'transparency'        0-1 transparency (0 is max transparency) 
%
% Output:
%       
%       Plot of spectral data with buckets drawn in the desired colors and
%       transparency.
%
% Usage: 
%         
%         [cinfo] = customColormap(pvals,'colors',[0.5686,0.0392,0.0392;1,1,1],'setPoints',[0,max(pvals)]);
%         colorBuckets(matrix,ppm,buckets,cinfo,'p-values');
%                 
% MTJ 2020


%%

%lineColors = lines(size(X,1));

%%
%     if ~isempty(varargin)   
%         ind = find(strcmp('lineColors',varargin),1);
%         if ~isempty(ind)
%             lineColors = varargin{ind+1}; % reset resolution to passed val
%         end
% 
%     end
        
%%
    % Several components as bars below would be nice
        %colors = [0.5686    0.0392    0.0392;1.0000    1.0000    1.0000;0  0  0.5686];
        
    % Generate the colormap
            
        %rng = max([min(abs(vals)),max(abs(vals))]);
        %[cinfo] = customColormap(vals,'setPoints',[-rng,0,rng]);

    % Plot the bins on the data
        
        figure, hold on
            highlightROIs(buckets',max(X(:)),'colorVect',bucketColors.rgb,'edgeColor','none','transparency',transparency)
            for i = 1:size(lineColors,1)
                plot(ppm,X(i,:),'Color',lineColors(i,:))
            end
            ax = gca;
            set(gca,'XDir','reverse')
            title({titleStr})
            xlabel('Chemical Shift (ppm)')
            ylabel('Signal Intensity')
            yticks(gca,[])
            addReasonableLegend(unique(lineLabels),unique(lineColors,'rows'))
            
            colormap(1 - (transparency * (1 - bucketColors.cmap))); % transparency applied mathematically
            t=colorbar;
            set(t,'Ticks',(bucketColors.setPts-min(bucketColors.setPts))/range(bucketColors.setPts))
            set(t,'TickLabels',num2cell(round(bucketColors.setPts,2)))
            set(get(t,'ylabel'),'String', colorbarStr)
            %set(t,'alpha',transparency)
            %alphamap(transparency);
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