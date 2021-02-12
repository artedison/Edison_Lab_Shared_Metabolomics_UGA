function [selectedRidgeInds,logicalInds] = clickRidge(ridges,rad,objType)

% Meant to be a more user-friendly version of selectLine().
% Picks closest line object to mouse click. 
% Call on active figure with line objects (objType = 'Line' is default,
% behavior for other graphics objects may be unpredictable)
% Press any key other than Enter/Return to click another ridge. 
% Mainly used in CIVM ridge tracking/tracing workflows.
% 

if ~exist('objType','var')
    objType = 'Line';
end

    % Initialize
    
    fig = gcf;
    
        % Keep track of selected ridges (index: ridge objects on plot)
            % Get line objects
            
                lineObjs = findall(gca,'Type',objType);
                
%                 rlengths = cellfun(@length, {ridges.ppms});
%                 
                % Get the unlisted inds
                
                    ridgeInds = {ridges.ppms};
                    ridgeInds = {lineObjs(:).XData};
                    for i = 1:length(ridgeInds)
                        ridgeInds{i} = ones(1,length([ridgeInds{i}]))*i;
                    end
                    
                    ridgeInds = [ridgeInds{:}]';
                    
                    lineData = [lineObjs(:).XData;lineObjs(:).YData]';
                    
%                 % Compare to make sure they're positive matches
%                 
                    ridgeCp = [{ridges.ppms};{ridges.times}];
                    id = false(length(lineObjs));
                    lineCp = [{lineObjs.XData};{lineObjs.YData}];
                    
                    ridgeCplens = cellfun(@length, {ridges.ppms});
                    lineCplens = cellfun(@length, {lineObjs.XData});
                    [jvec,ivec] = find(ridgeCplens == lineCplens');
                    
                    for i = 1:length(ivec)
                            eq = [lineCp{:,jvec(i)}] == [ridgeCp{1,ivec(i)},ridgeCp{2,ivec(i)}']';
                            id(ivec(i),jvec(i)) = all(eq(1:(size(eq,1)+1):end)); % always has a diagonal
                    end
                    [~,LmatchR] = find(id);
                    lineObjs = lineObjs(LmatchR);
                    
        % Keep track of how many times a ridge was clicked (odd = selected)
        
            % start as zeros
            
                bold = zeros(length(lineData),1);

     % Look for rad (radius = maximum distance)
    
        if ~exist('rad','var')
            d = pdist2(lineData,lineData,'Euclidean');
            rad = max(d(:));
        end               
                
    % While ~enter
        select = 1;
    
    
    while select
        
        % Prompt for click or button press
        
            figure(fig)
            title('Press any letter key to select another ridge, or press Enter/Return to exit')
            
            waitforbuttonpress;
                key = get(gcf,'CurrentCharacter');
                
                if ~(key == 13)
                
                    % Click mouse on points on the plot
                    
                        [x,y] = ginput(1);
                        
                    % Within a given radius (euclidean distance)
                        
                        % Which ridge has points closest to the cursor?
                            
                            distances = sqrt( sum(([x,y] - lineData) .^ 2, 2) );
                            [val,closestPtInd] = min( distances );
                            
                            % Is any ridge selected within radius?
                            
                            if val<=rad
                                selectedRidge = ridgeInds(closestPtInd);
                            else
                                % if not, then jump to next iteration
                                selectedRidge = [];
                                continue
                            end

                        % If zero, flip to one. If 1, flip to zero
                            
                            bold(selectedRidge) = abs(1-bold(selectedRidge));

                        % If 1, update object to be bolded
%                         for i = 1:length(lineObjs)
%                             lineObjs(i).LineWidth = 5;
%                         end
                        
                            if bold(selectedRidge)
                                lineObjs(selectedRidge).LineWidth = 9;
                            else
                                lineObjs(selectedRidge).LineWidth = 5;
                            end
                        
                else
                    % If buttonpress was 'Enter/Return', break out of the
                    % loop
                    break
                end
    end
    
    % Report the ridge inds

        logicalInds = bold;
        selectedRidgeInds = find(bold);
            
end