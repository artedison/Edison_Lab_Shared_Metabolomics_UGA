function [uniq_colors,c_per_row]= plotr_colorby(ppm,matrix,metadata,varargin)

%% Goncalo Gouveia Dec-5-2019
%function to plot and color acording to a metadata vector (color by
%group/replicate/etc)
%Also creates a legend with equivalent colors to the spectra

%input
%ppm - a ppm vector from NMR spectra
%matrix - any matrix where rows are samples and columns are matched to the
%ppm vector
%metadata - a vector of same size as the number of rows in matrix - can be
%a number or string - note that the legend function is capped at 50
%entries, more than that and they will not be displayed 
% [uniq_colors,c_per_row] outputs - uniq_colors, gives a matrix of three
% colums (each for an rgb value, where the rows correspond to the unique
% number of groups. c_per_row - outputs a matrix of 3 columns (each an rgb
% value, and the rows correspond to all the samples (samples of the same
% group have the same color)

%% Goncalo Gouveia Jan-7-2020
%updated to produce the color vectors as outputs (GG and MJ)

%Uses functions: 
    %Matlab version: 2019a (older versions may not work)
    %Matlab: string - ismember - numel - plot
    %Edisonlab functions: distinguishable_colors - plotr
    

ncolrs = string(unique (metadata,'rows')); %determine the number of unique elements in the metadata and transform into a string 

uniq_colors = distinguishable_colors(numel(ncolrs)); % creates a color map for the number of elements in the metadata

[LIA,LOCB]  = ismember (string(metadata),ncolrs,'rows'); %compares the metadata and the unique elements and get indices

c_per_row= uniq_colors(LOCB,:);%mapps colors to indices 



co = get(gca,'ColorOrder'); % Initial;
% Change to new colors.
set(gca, 'ColorOrder', c_per_row, 'NextPlot', 'replacechildren');
co = get(gca,'ColorOrder'); % Verify it changed;
% Now plot with changed colors.
hold on
plotr (ppm,matrix,varargin{:})

hold on

legend off
clear L

for i=1:size(uniq_colors,1)
 hold on   
L (i)= plot ( [NaN,NaN], 'color', uniq_colors(i,:), 'DisplayName',string(ncolrs(i)) );%creates a fake plot of nan data to get the unique values/colors for the plot

hold off
end

legend (L,(ncolrs),'NumColumns',2,'Box','off')
hold off
hold off
highlightLine(string(metadata))
hold on

% 