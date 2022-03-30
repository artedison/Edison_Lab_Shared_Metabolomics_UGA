function [ROIs,patches,fig] = extractROIs(varargin)
%% extractROIs

% Author: MTJ
% Version: 0.1
% Date: 2020
%
% Description:
%
%     Called with an open or saved figure in the working directory, and
%     extracts any buckets that are drawn as patch objects on that figure
%     to be used in other processes. Like reading a mask from the figure. 
%     Vectorized (faster). ROI(s) = Region(s) Of Interest.
% 
% Inputs:
%
%     Optional : 
% 
%         figure name or figure handle (if accessing saved figure without
%         opening)
%
% Output:
%       
%      ROIs     bucket/ROI boundaries 
%      patches  handles to patch objects in figure (in case modification is
%               desired)
%      fig      figure handle
%
% Usage: 
%         
%         [ROIs,patches,fig] = extractROIs();
%         [ROIs,patches,fig] = extractROIs('My Bucket Figure.fig');
%         [ROIs,patches,fig] = extractROIs(fig);
%                 
% MTJ and GG 2020


    if ~isempty(varargin)
        
        % Check to see if arg was a figure name
            type = cellfun(@class,varargin,'UniformOutput',false);
            isname = strcmp(type,'char');
            if any(isname)
                open(varargin{find(isname,1)});
                fig = gcf;
            end
            
        % Check to see if arg was a figure handle (supersedes fig name)
            isfig = strcmp(type,'matlab.ui.Figure');
            if any(isfig)
                fig = varargin{find(isfig,1)};
            end
        
    else
        fig = gcf; % default is to get the current figure
    end
    
        D=get(gcf,'Children'); %get the handle of the line object
        patches = findall(D,'Type','Patch');  % make a list of all patch boxes in the figure
        if ~isempty(patches)
            verts = [patches.Vertices];
            ROIs = flipud([verts(3:8:(numel(verts)));verts(1:8:(numel(verts)))])';
        else
            ROIs = [];
        end
end