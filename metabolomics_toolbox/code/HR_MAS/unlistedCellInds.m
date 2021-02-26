function [contents,inds,inds_cells] = unlistedCellInds(cells)

% Convert cells of double arrays to cells of double
% arrays containing the indices of the cells:
% E.g. 
% cells = [ {[1,2,3]}, {1,2,6,4,3]}, {3,2,1,9,5}];
%   -> contents:    [1,2,3,1,2,6,4,3,3,2,1,9,5];
%   -> inds:        [1,1,1,2,2,2,2,2,3,3,3,3,3];
%   -> inds_cells:  [ {[1,1,1]}, {2,2,2,2,2]}, {3,3,3,3,3}];
%
% MTJ 2021

        inds_cells = cells;
        
        for i = 1:length(inds_cells)
            inds_cells{i} = ones(1,length([inds_cells{i}]))*i;
        end

        inds = [inds_cells{:}]';

        contents = [cells{:}]';

end