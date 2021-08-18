function [lineInfo,lineObjs] = extractLineData()

% Extract line data (e.g. for ridge selection)
% Simplifies vectorized operations on ridges. 
% Only works in 2D right now. Could add zdata pretty easily, though.
% no inputs
% outputs:
%   lineInfo.unlisted.lineInds: line number for each row (point) in lineInfo.unlisted.lineData
%   lineInfo.unlisted.lineData: (:,1) is x data, (:,2) is y data for the lines
%   lineInfo.XData:   x data in form of cells (corresponding to lines)
%   lineInfo.YData:   y data in form of cells (corresponding to lines)
%
%
% MTJ 12FEB2020

    % Get the line data
    
            lineObjs = findall(gca,'Type','Line');
            if isempty(lineObjs)
                error(['extractLineData: It looks like the figure this was called on doesn''t contain objects of type ''Line''. Usually, this means there was no figure open when the function was called.'])
            end

        % Get the unlisted inds for line obj points (after
        % reordering, to compare with mouse clicks)

            lineInds = {lineObjs(:).XData};
            
            for i = 1:length(lineInds)
                lineInds{i} = ones(1,length([lineInds{i}]))*i;
            end

        % Unlist the data
            lineInfo.unlisted.lineInds = [lineInds{:}]';
            lineInfo.unlisted.lineData = [lineObjs(:).XData;lineObjs(:).YData]';

        lineInfo.XData = lineObjs.XData;
        lineInfo.YData = lineObjs.YData;


end