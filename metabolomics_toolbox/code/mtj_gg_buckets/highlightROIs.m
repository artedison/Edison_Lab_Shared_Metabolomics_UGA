function [] = highlightROIs(ROIs,height,varargin)

%% highlightROIs

% Author: MTJ
% Version: 0.3
% Date: 2020
%
% Description:
%
%   'Highlights' the regions of a plot provided in 'ROIs'
%   (typically ppm values on a spectral plot, for plotting bins/buckets).
%   Draws the regions as patch objects of specified height and boundaries
%   provided in ROIs. Draws on active plot, or generates a new figure and
%   plot. Transparency creates the highlighting effect. Shading is additive
%   in regions of overlap between ROIs. Functionality is provided for
%   highlighting ROIs on stackplots (e.g. made with the stackSpectra()
%   function.
%
% Inputs:
%
%       Required Arguments:
%             ROIs:         Regions Of Interest. 2 x n matrix of x axis values
%                           (e.g. ppm values), where the ROIs(1,n) is the lower
%                           bound and ROIs(2,n) is the upper bound for ROI n. 
% 
%             height:       height of the highlights/patch objects. Uniform for 
%                           all patches.
%
%                   All of the following are currently required for 
%                   stackplots (e.g. stackSpectra() plots), but are not
%                   used in other cases:
% 
%                     'horzshift'           horizontal shift factor used in 
%                                           call to stackSpectra()
%                     'vertshift'           vertical shift factor used in 
%                                           call to stackSpectra()
%                     'numberOfSpectra'     number of spectra in the 
%                     'extension'           factor controlling how far the
%                                           highlight extends on each end, 
%                                           since highlights on slanted 
%                                           spectra appear shorter.
%
%       Optional Name,Value pair Arguments: 
%           
%             'color'               base highlight color (not considering
%                                   transparency). RGB or MATLAB string
%                                   definition (e.g. 'k'). Default 'r'.
%             'transparency'        degree of highlight transparency. 
%                                   Passthrough to the 'FaceAlpha' 
%                                   property of the Patch objects. Values 
%                                   0 to 1. Lower is more transparent. 
%             'edgeColor'           color of the border for the regions.
%                                   Defined as in 'color'. Default 'none'.
% 
%
% Output:
%       
%       Patch objects corresponding to ROIs on the current active axis.
%       These will be drawn on top of existing objects in the plot.
%
% Usage: 
%         
%       highlightROIs(ROIs,max(matrix(:)),varargin)
%       highlightROIs(ROIs,max(matrix(:)),'edgeColor','k','color',[0.5,0.75,0.3],'transparency',0.05)
%                 
% MTJ 2017

%%
edgeColor = 'none';

    % Read in any options:
        if ~isempty(varargin)
            varnames = varargin(1:2:end);
            varvals = varargin(2:2:end);

            varInd = contains(varnames,'color');
            if any(varInd)
                color = varvals{varInd};
            end
            varInd = contains(varnames,'transparency');
            if any(varInd)
                transparency = varvals{varInd};
            end            
            varInd = contains(varnames,'edgeColor');            
            if any(varInd)
                edgeColor = varvals{varInd};
            else
                edgeColor = 'none';
            end                
            
            
            % StackSpectra params:
            varInd = contains(varnames,'horzshift');
            if any(varInd)
                stackParams.horzshift = varvals{varInd};
            end
            varInd = contains(varnames,'vertshift');
            if any(varInd)
                stackParams.vertshift = varvals{varInd};
            end    
            varInd = contains(varnames,'numberOfSpectra');
            if any(varInd)
                stackParams.numberOfSpectra = varvals{varInd};
            end    
            varInd = contains(varnames,'extension');
            if any(varInd)
                stackParams.numberOfSpectra = varvals{varInd};
            end    
        end
        
    % Make the highlights
        for i = 1:size(ROIs,2)
            % Calculate the coordinates
                width = (ROIs(2,i)-ROIs(1,i));
                llx = ROIs(1,i); % 
                lux = llx;
                rlx = llx + width; % = rux
                rux = llx + width; % = rux;
                lly = 0;
                luy = height;
                rly = 0;
                ruy = height;
                xcoords = [llx rlx rux lux];
                ycoords = [lly rly ruy luy];
                
            % If it's a stackSpectra plot, modify the coordinates:
                if exist('stackParams','var')
                    
                    % luy and ruy coords need to shift based on vertshift, numberOfSpectra
                        maxVertShift = stackParams.vertshift * stackParams.numberOfSpectra;
                        luy = 0 + stackParams.vertshift * stackParams.extension;
                        ruy = 0 + stackParams.vertshift * stackParams.extension;
                        lly = lly - maxVertShift - stackParams.vertshift * stackParams.extension;
                        rly = rly - maxVertShift - stackParams.vertshift * stackParams.extension;
                        ycoords = [lly rly ruy luy]; 
                        
                    % lux and rux coords need to shift based on horzshift, numberOfSpectra
                        maxHorzShift = stackParams.horzshift * stackParams.numberOfSpectra;
                        llx = llx + stackParams.horzshift * stackParams.extension;
                        rlx = rlx + stackParams.horzshift * stackParams.extension;
                        lux = lux - maxHorzShift;
                        rux = rux - maxHorzShift;
                        xcoords = [llx rlx rux lux]; 
                end

            % Make the patch:
                p=patch(xcoords,ycoords,'r');

            % Modify the patch in various ways:
                % Face Color
                    if exist('color','var')
                        %p=patch(xcoords,ycoords,color); %light red % pink [1    0.5  1]
                        set(p,'FaceColor',color);
                    else
                        %p=patch(xcoords,ycoords,'r'); %light red % pink [1    0.5  1]
                        set(p,'FaceColor','r'); %light red % pink [1    0.5  1]
                    end

                % Face Color
                    if exist('edgeColor','var')
                        %p=patch(xcoords,ycoords,color); %light red % pink [1    0.5  1]
                        set(p,'EdgeColor',edgeColor);
                    else
                        %p=patch(xcoords,ycoords,'r'); %light red % pink [1    0.5  1]
                        set(p,'EdgeColor','none'); %light red % pink [1    0.5  1]
                    end
                    
                % Transparency
                    if exist('transparency','var')
                        set(p,'FaceAlpha',transparency)
                    else
                        set(p,'FaceAlpha',0.1);
                    end

        end
end