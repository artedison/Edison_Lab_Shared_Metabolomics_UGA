function [selectedRidgeInds,logicalInds,ridgeInds,lineObjs,lineData] = clickRidge(ridges,titleStr,rad,objType)

% Meant to be a more user-friendly version of selectLine() that also avoids
% global variables (often a poor coding practice).
% Picks closest line object to mouse click. 
% 
% Call on active figure with line objects (objType = 'Line' is default,
% behavior for other graphics objects may be unpredictable)
% Press any key other than Enter/Return to click another ridge. 
% Mainly used in CIVM ridge tracking/tracing workflows.
% Ridgepoints will be sorted by time before being compared for exact match.
%
% Params:
% 'ridges' struct array must contain the following for each cluster i:
%     ridges(i).ppms = [];
%     ridges(i).times = [];
%
%   this is required to match the list of objects in the figure with an
%   the list you have in-hand. 
%
% Optional:
% rad is the max search distance (set automatically, but can adjust)
% objType allows different objects (e.g. scatter points) to be clicked. 
%     'Line' is default. Currently not tested. 
% 
% 
% MTJ 12FEB2021
% 

    % Initialize
    
        % Maybe sometime later, this can be extended/adapted to diverse objects
        
            if ~exist('objType','var')
                objType = 'Line';
            end

        % In case another title should be passed
            
            if ~exist('titleStr','var')
                titleStr = '';
            end
        
        % Get the figure handle to make sure it can be called to the top later
        
            fig = gcf;
    
        % Get line objects

                lineObjs = findall(gca,'Type',objType);
                if isempty(lineObjs)
                    error(['clickRidge: It looks like the figure this was called on doesn''t contain objects of type ''',objType,''', so we can''t select them. Usually, this means there was no figure open when the function was called.'])
                end
                                                     
       % Make sure line objs match ridges. It is assumed that no
       % lines exist without ridges, but not vice-versa.

            % We're interested in comparing the ppm and time vals
            % for the ridges. (ppm,time) maps 1to1 to surface intensity, so intensity is unnecessary

                ridgeCp = [{ridges.ppms};{ridges.times}];
                    % Sort each ridge's points by time
                    for i = 1:size(ridgeCp,2)
                        [ridgeCp{2,i},si] = sort(ridgeCp{2,i});
                        ridgeCp{1,i} = ridgeCp{1,i}(si);
                    end
                lineCp = [{lineObjs.XData};{lineObjs.YData}]; 
                    % Sort each line's points by time
                    for i = 1:size(lineCp,2)
                        [lineCp{2,i},si] = sort(lineCp{2,i});
                        lineCp{1,i} = lineCp{1,i}(si);
                    end
                    
            % Need an empty matrix to hold all pairwise comparisons between cell contents from the two lists of cells (concat ppms and times)
                id = false(length(lineObjs), length(ridges)); 

            % Find all pairwise length matches (to eliminate most unnecessary comparisons)
                ridgeLengths = cellfun(@length, {ridges.ppms});
                lineLengths = cellfun(@length, {lineObjs.XData});
                [jvec,ivec] = find(ridgeLengths == lineLengths'); % use these inds to drive a for loop

                for i = 1:length(ivec)
                    
                     % Compare concatenation of x positions and y positions (z is implied, unnecessary)
                        eq = [lineCp{:,jvec(i)}] == [ridgeCp{1,ivec(i)},ridgeCp{2,ivec(i)}']';

                     % If the ridges are actually the same, then the diagonal will be all ones
                        id(ivec(i),jvec(i)) = all(eq(1:(size(eq,1)+1):end)); % always square matrix bc we only compare vectors of matching lengths (filter above)
                end
                    
        % Reorder the line objects within the function to match
        % the provided ridges. Two reasons to do this:
        %     1) filter out any line objects without ridges
        %     2) set lineObjs handle order to match ridge order

            [RmatchL,~] = find(id);         % line objects are on the columns
            lineObjs = lineObjs(RmatchL);   
                    
                % Get the unlisted inds for line obj points (after
                % reordering, to compare with mouse clicks)
                
                    ridgeInds = {lineObjs(:).XData};
                    for i = 1:length(ridgeInds)
                        ridgeInds{i} = ones(1,length([ridgeInds{i}]))*i;
                    end
                    
                    ridgeInds = [ridgeInds{:}]';
                    
                    lineData = [lineObjs(:).XData;lineObjs(:).YData]';
                    
                    
        % Keep track of how many times a ridge was clicked (odd = selected)
        
            % start as zeros
            
                bold = zeros(length(lineObjs),1);

     % Look for rad (radius = maximum distance)
    
        if ~exist('rad','var')
            d = pdist2(lineData,lineData,'Euclidean');
            rad = max(d(:));
        end               
                
    % While ~enter
        select = 1;
        
%% New method of ridge selection 
% via https://www.mathworks.com/matlabcentral/answers/325754-how-to-stop-while-loop-by-right-mouse-button-click
% This way, we just click until a button press kicks us out of the loop
%

    ax = gca;
    fig = ancestor(ax, 'figure');
    hold(ax, 'on');
    
    while true
        
            % Click mouse on a point on the plot

                [x,y] = ginput(1);
                
            % If buttonpress was 'Enter/Return', break out of the
            % loop

            % Figure out if the action was a click or a buttonpress
            
                sel = get(fig, 'SelectionType');

                if strcmpi(sel, 'alt')
                    break
                end

            % Within a given radius (euclidean distance)

                % Which ridge point is closest to the cursor?

                    distances = sqrt( sum(([x,y] - lineData) .^ 2, 2) );
                    [val,closestPtInd] = min( distances );

                    % Is any ridge selected within radius?

                        if val<=rad
                            % If so, which ridge does it belong to?
                                selectedRidge = ridgeInds(closestPtInd);

                            % If zero, flip to one. If 1, flip to zero

                                bold(selectedRidge) = abs(1-bold(selectedRidge));

                            % If 1, update object to be bolded

                                if bold(selectedRidge)
                                    lineObjs(selectedRidge).LineWidth = 9;
                                else
                                    lineObjs(selectedRidge).LineWidth = 5;
                                end

                        else
                            % if not, then jump to next iteration
                            % (selectedridge stays set to [])
                            continue
                        end

    end
    
    hold(ax, 'off')   

    
%% Requires button press between each ridge and ends using "Enter/Return"
%
%     while select
%         
%         % Prompt for click or button press
%         
%             figure(fig)
%             title({titleStr,'Press any letter key to (de-)select another ridge, or press Enter/Return to exit'})
%             
%             waitforbuttonpress;
%                 key = get(gcf,'CurrentCharacter');
%                 
%                 if ~(key == 13)
%                 
%                     % Click mouse on points on the plot
%                     
%                         [x,y] = ginput(1);
%                         
%                     % Within a given radius (euclidean distance)
%                         
%                         % Which ridge point is closest to the cursor?
%                             
%                             distances = sqrt( sum(([x,y] - lineData) .^ 2, 2) );
%                             [val,closestPtInd] = min( distances );
%                             
%                             % Is any ridge selected within radius?
%                             
%                             if val<=rad
%                                 % If so, which ridge does it belong to?
%                                     selectedRidge = ridgeInds(closestPtInd);
% 
%                                 % If zero, flip to one. If 1, flip to zero
% 
%                                     bold(selectedRidge) = abs(1-bold(selectedRidge));
% 
%                                 % If 1, update object to be bolded
%                         
%                                     if bold(selectedRidge)
%                                         lineObjs(selectedRidge).LineWidth = 9;
%                                     else
%                                         lineObjs(selectedRidge).LineWidth = 5;
%                                     end
%                             
%                             else
%                                 % if not, then jump to next iteration
%                                 % (selectedridge stays set to [])
%                                 continue
%                             end
%                         
%                 else
%                     % If buttonpress was 'Enter/Return', break out of the
%                     % loop
%                     break
%                 end
%     end
%%    
    % Report the ridge inds

        logicalInds = bold;
        selectedRidgeInds = find(bold);
            
end