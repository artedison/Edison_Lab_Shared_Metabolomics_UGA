function addReasonableLegend(labels,colors,varargin)
%% addReasonableLegend

% Author: MTJ, GG
% Version: 0.1
% Date: 2020
%
% Description:
%
%   Add single legend for each sample type (label) instead of for each
%   individual sample. Stolen from GG plotr_colorby and functionalized for
%   use in other contexts. GG, you're a legend for providing this solution.
%
% Inputs:
%
%     Required Arguments:
%         labels:       cell array or numeric array defining n sample group identities. These will be the legend entries.
%         colors:       nx3 numberic array defining rgb colormap (each row contains [r,g,b] values [0,1] for each label, respectively)
%
%     Optional Name,Value pair Arguments: 
% 
%         'addBox'      any value works; if set to anything, then a box is drawn around the legend
%         'numColumns'  must be a whole number. Number of legend columns you want
%         'textSize'    number defining text size in legend
%         'lineWidths'  width of the lines (to match your spectra)
%
% Output:
%       A legend on the current active plot corresponding to the categories provided.
%       
%
% Usage: 
%         addReasonableLegend(labels,colors)
%         addReasonableLegend(labels,colors,'addBox','true')
%         addReasonableLegend(labels,colors,'numColumns',2)
%         addReasonableLegend(labels,colors,'addBox','true','numColumns',2,'textSize',15,'lineWidths',[2,1])
%                 
% MTJ and GG 2020

%% Parse Optional Params

% Following name/value pair format. 
if ~isempty(varargin)
    for k = 1:2:length(varargin)
        switch varargin{k}
            
            case 'addBox'
                addBox = varargin{k+1};
                
            case 'numColumns'
                numColumns = varargin{k+1};
                % Check to make sure numColumns is appropriate
                if ~isnumeric(numColumns)
                    error('ERROR: addReasonableLegend: numColumns argument must be a whole number > 0.')
                end
                
            case 'textSize'
                textSize = varargin{k+1};
                
            case 'lineWidths'
                lineWidths = varargin{k+1};
                if length(lineWidths)<length(labels)
                    lineWidths = repmat(lineWidths,1,length(labels));
                end
        end
    end
end

% Set default values
if ~exist('numColumns','var')
    numColumns = 1;
end
if ~exist('textSize','var')
    textSize = 15; 
end
if ~exist('lineWidths','var')
    defaultLineWidth = 1;
    lineWidths = repmat(defaultLineWidth,1,length(labels)); 
end

%% Make the legend

    %legend off
    %figure(axesHandle)
    
    hold on
    for i=1:size(colors,1)
        L(i)= plot ( [NaN,NaN], 'color', colors(i,:), 'DisplayName',string(labels(i)),'LineWidth',lineWidths(i));%creates a fake plot of nan data to get the unique values/colors for the plot
    end
    
    if exist('addBox','var')
        lgd = legend (L,(labels),'NumColumns',numColumns,'Box','on');
    else
        lgd = legend (L,(labels),'NumColumns',numColumns,'Box','off');
    end
    
    lgd.FontSize = textSize;
    set(lgd, 'Interpreter', 'none')
    
end