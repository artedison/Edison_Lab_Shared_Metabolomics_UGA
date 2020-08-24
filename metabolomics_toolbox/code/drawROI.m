function [xbds,ybds] = drawROI() 
%% drawROI
%     Draw a rectangle on a figure to get the x and y boundaries. Useful for 
%     defining a Region Of Interest (ROI) for end removal, solvent peak 
%     removal, or integration region definition.
%     
%     Inputs: none
%     
%     Outputs: 
%             xbds    x coordinates of ROI (lowest to highest)
%             
%             ybds    y coordinates of ROI (lowest to highest)
% 
%             
%
% MTJ 2020

        rect = getrect(gca);
            xpos = [rect(1),rect(1)+rect(3)];
            ypos = [rect(2),rect(4)];
            xbds = xpos;
            ybds = ypos;
            % Sort ppms (with intensities) of chosen point pair so
            if xbds(2) < xbds(1)
                xbds = flip(xbds);
            end
            if ybds(2) < ybds(1)
                ybds = flip(ybds);
            end
end