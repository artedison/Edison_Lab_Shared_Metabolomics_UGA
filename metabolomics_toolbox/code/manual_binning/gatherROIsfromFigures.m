function [features] = gatherROIsfromFigures(listOfFigures,matrix,ppm,peakpicking,method,shifting)
%% Function for gathering ROIs as pink patch boxes from a list of figures
%{
    listOfFigures: cell array of figure names (must be in the working
    directory)

    matrix:     the spectral matrix provided for manual binning; the one in
                    the figures
    ppm:        the ppm vector ""
    peakpicking: option to pick a peak in each ROI. Useful for reference to
                    each region (default); '' is default
                    'off' disables peakpicking
    method:     peakpicking methods for basicPeakPickingROIs.m:
                    'max' picks the highest point in each region (default)
                    'middle' picks the middle of the region as the nominal
                        "peak"                    
%}
%% Import Regions from patch boxes in saved figures

    oldfigs = listOfFigures; % list of figure names saved above
    oldboxes = []; % initializing

    set(0,'DefaultFigureVisible','off'); % turn figures off to save time
        for j = 1:length(oldfigs) % loop through the figures in the list
            f = open(oldfigs{j}); 
            fprintf(['\nGetting ROIs from "',oldfigs{j},'"...\n']);
            %% Get the current figure's data
                %open data.fig %open your fig file, data is the name I gave to my file
                D=get(gcf,'Children'); %get the handle of the line object
                oldRecs = findall(D,'Type','Patch');  % make a list of all patch boxes in the figure
                    for i = 1:length(oldRecs)       % for each patch box, get the vertices and store
                        verts = oldRecs(i).Vertices;
                        oldboxes = horzcat(oldboxes,[verts(3);verts(1)]);
                    %{
                        verts = [1 2; % ll % Only need one l and one r; 1&3
                                 3 4; % lr
                                 5 6; % ur
                                 7 8];% ul
                    %}
                    end
            close(f)
        end
        
            set(0,'DefaultFigureVisible','on'); % turn figures back on
            
        if size(oldboxes,2)>0

               % Sort features by ppm value from right (small) to left (large)
                        %oldboxes = [oldboxes,y]; % from cleanup run
                    % Make sure smaller ppm is on top:
                        oldboxes = sort(oldboxes,1);
                    % Sort features high to low; store in features_manual:
                        oldboxes(3,:) = oldboxes(1,:) + oldboxes(2,:); % get sum of the ppms
                        oldboxes = sortrows(oldboxes',3);              % sort by sum, the idea being that
                                                                        % if two features overlap, the one
                                                                        % that reaches farther will be second.
                        features_manual = oldboxes(:,1:2)';             % store it


            %% Bin each region

                % Initialize features structure
                    features = struct();
                    %features.leftBounds = features_manual(2,:);
                    %features.rightBounds = features_manual(1,:);
                    %features.bounds = features_manual;
                    features.bounds = ppm(matchPPMs(features_manual,ppm));
                    features.intensities = nan(size(matrix,1),size(features.bounds,2));
                    [indices_manual] = matchPPMs(features_manual,ppm);      % convert feature bounds (ppm) to ppm indices
                % Integrate each region and store in features structure alongside
                % bounds
                    for i = 1:size(matrix,1)
                        for j = 1:size(indices_manual,2)
                            features.intensities(i,j) = sum(    matrix(i,   indices_manual(1,j):indices_manual(2,j))    , 2);
                        end
                    end

                %% Get a list of ppm indices for each feature
                        % concatenate for use as a single list of indices for use in eval(sprintf('subsetColumns = matrix(:,%s);',features.regionsIndicesList))
                    regionsIndices = [num2str(matchPPMs(features.bounds(1,:)',ppm)),repmat(':',size(features.bounds,2),1),num2str(matchPPMs(features.bounds(2,:)',ppm))];
                    features.regionsIndicesList = ['[',regionsIndices(1,:)];
                        for i = 2:(size(regionsIndices,1)-1)
                            features.regionsIndicesList = [features.regionsIndicesList,',',regionsIndices(i,:)];
                        end
                    features.regionsIndicesList = [features.regionsIndicesList(find(~isspace(features.regionsIndicesList))),']'];

              if ~strcmp(peakpicking,'off')
                %% Peak pick highest peak in each ROI
                    %method = 'middle';
                    %method = 'max';
                    if ~(or(strcmp(method,'max'),strcmp(method,'middle')))
                        features.peaks = basicPeakPickingROIs(matrix,ppm,features,method,shifting);
                    else
                        features.peaks = basicPeakPickingROIs(matrix,ppm,features,'max',shifting);
                    end
                        fprintf(['\n\tGathered, binned, and peak-picked ',num2str(size(features.intensities,2)),' ROIs from ',num2str(length(oldfigs)),' figures\n\n']);
              else
                    fprintf(['\n\tGathered, and binned ',num2str(size(features.intensities,2)),' ROIs from ',num2str(length(oldfigs)),' figures\n\n']);

              end
%%                    

        else

            fprintf('\nReturning an empty features structure\n')
            %close(gcf)
            features = struct();
            features.bounds = [];
            features.intensities = [];
            features.regionsIndicesList = [];
            features.peaks = [];
        end
end