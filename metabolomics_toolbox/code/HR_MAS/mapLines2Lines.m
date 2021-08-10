function [lineObjs,lineObjData,inds_myLines,indsCells_myLines] = mapLines2Lines(lineObjs,myLines_x,myLines_y)

% figureLines = matlab figure line object array, i.e. from findall(gca,'Type','Line')
% myLines = 
% 
% MTJ 12FEB2021

                                                     
       % Make sure line objs match ridges. It is assumed that no
       % lines exist without ridges, but not vice-versa.

            % We're interested in comparing the ppm and time vals
            % for the ridges. (ppm,time) maps 1to1 to surface intensity, so intensity is unnecessary

                myLinesCp = [myLines_x;myLines_y];
                lineCp = [{lineObjs.XData};{lineObjs.YData}];  

            % Need an empty matrix to hold all pairwise comparisons between cell contents from the two lists of cells (concat ppms and times)
                id = false(length(lineObjs), length(myLines_x)); 

            % Find all pairwise length matches (to eliminate most unnecessary comparisons)
                ridgeLengths = cellfun(@length, myLines_x);
                lineLengths = cellfun(@length, {lineObjs.XData});
                [jvec,ivec] = find(ridgeLengths == lineLengths'); % use these inds to drive a for loop

                for i = 1:length(ivec)
                    
                     % Compare concatenation of x positions and y positions (z is implied, unnecessary)
                        eq = [lineCp{:,jvec(i)}] == [myLinesCp{1,ivec(i)},myLinesCp{2,ivec(i)}']';

                     % If the ridges are actually the same, then the diagonal will be all ones
                        id(ivec(i),jvec(i)) = all(eq(1:(size(eq,1)+1):end)); % always square matrix bc we only compare vectors of matching lengths (filter above)
                end
                    
        % Reorder the line objects within the function to match
        % the provided ridges. Two reasons to do this:
        %     1) filter out any line objects without ridges
        %     2) set figureLines handle order to match ridge order

            [RmatchL,~] = find(id);         % line objects are on the columns
            lineObjs = lineObjs(RmatchL);   
                    
                % Get the unlisted inds for line obj points (after
                % reordering, to compare with mouse clicks)
                
                    [~,inds_myLines,indsCells_myLines] = unlistedCellInds({lineObjs(:).XData});
                    lineObjData = [lineObjs(:).XData;lineObjs(:).YData]';

end