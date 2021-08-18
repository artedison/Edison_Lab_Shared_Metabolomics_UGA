function [xbds,ybds] = selectRegion()
%% selectRegion
% Select a region from a plot and return the coordinates
% 
%     Inputs:     none, but have an active figure open
% 
%     Outputs:    
%             xbds    x coordinates of ROI (lowest to highest)
%             
%             ybds    y coordinates of ROI (lowest to highest)
% 
%
% MTJ 2020

    % Make a button to activate the selection (e.g. after pan/zooming)
    
        uicontrol('Position',[2 4 300 30],'String','Click to draw rectangle around ROI.',...
                  'Callback','uiresume(gcbf)');

        uiwait(gcf);% Wait for input, continue code when button is clicked
    
    % Select
    
        [xbds,ybds] = drawROI();

end