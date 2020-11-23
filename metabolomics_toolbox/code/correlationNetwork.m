function [Rmat,Pmat,RmatT,PmatT,interactionsTable] = correlationNetwork(matrix,Rthresh,Pthresh,filename,labels,RmatIn,PmatIn,removeSelfEdges,whichCorr,additionalColumns)
%% correlationNetwork

%   Author: MTJ
%   Version: 0.1
%   Date: 2017
%
%   Description 
%
%       Takes or creates a correlation matrix (optional), applies hard cutoff
%        (threshold) for correlation value or p value, then builds an
%        interaction table for visualization in network visualization software 
%        like Cytoscape. Mainly used as a formatting function for that file, and
%        also applies common filters (e.g. remove duplicate edges, remove self
%        edges). 
% 
%   Inputs
%
%         matrix:           dataset where feature (variable/node) number corresponds to 
%                           column and sample corresponds to row.
%         Rthresh:          the Pearson Correlation Coefficient cutoff (everything below is
%                           set to R = 0
%         Pthresh:          cutoff for the PCC p-value, as given by corr() in Matlab. All
%                           p-values greater than the cutoff are ignored.
%                   NOTE: The intersection of edges passing Rthresh and Pthresh is given 
%                       as a result. 
%                   NOTE: The pvals for the diagonal are always 1, but we change this
%                       to zero in the rare event that we want to preserve
%                       self-correlations. For reference, see Matlab corr(), lines 317-328.
%         filename:         output filename (as a string, .csv included). This is placed
%                           wherever the current wd is.
%         labels:           option to pass specific column (feature) labels as cell
%                           array of strings or a double/int array. I usually just use 
%                           the feature index (double array 1:size(matrix,2)).
%         RmatIn:           typically []; this is for the case where you want to bypass
%                           the distance calculation (e.g. use the function for its 
%                           filters/formatting/file output). Contains the correlation
%                           values (e.g. Pearson's Rho) in an (m x m) matrix where the
%                           initial data had n measurements of m variables. 
%         PmatIn:           typically []; the pvals corresponding to RmatIn.
%         removeSelfEdges:  removes self-edges if set to 'removeSelfEdges',
%                           otherwise autocorrelations are preserved.
%                    NOTE: The pvals for the diagonal are always 1, but we change this
%                       to zero in the rare event that we want to preserve
%                       self-correlations. For reference, see Matlab corr(), lines 317-328.
%         whichCorr:        allows for filtering for positive correlations only
%                           'positive', negative correlations only 'negative', or
%                           both (default, 'both')
%   Outputs
%
%         Rmat:                 correlation matrix
%         Pmat:                 p value matrix
%         RmatT:                thresholded correlation matrix (R < Rthresh = 0)
%         PmatT:                thresholded p value matrix (p > Pthresh = 1)
%         interactionsTable:    edge table (node, node, correlation)
%
% 
%
%   Usage:
%             filename = 'TestFeatures_Network_AlltoAll.csv';
%             Rthresh = 0.9;
%             Pthresh = 0.05;
%             matrix = features.intensities; % typically binned data
%             labels = 1:size(matrix,2); % default is feature number
%             RmatIn = Rmat;%[];
%             PmatIn = Pmat;%[];
%
%             
%             [Rmat,Pmat,RmatT,PmatT,interactionsTable] = correlationNetwork(matrix,Rthresh,Pthresh,filename,labels,RmatIn,PmatIn,'removeSelfEdges','both');
%

tic
%% Do the correlation calculations:

    if or(isempty(RmatIn),isempty(PmatIn))
        fprintf('\nCalculating correlations...\n')
        [Rmat,Pmat] = corr(matrix); % for a matrix (m x n), Rmat and Pmat are both (m x m).
                                    % *** pvals on the diagonal, while
                                    % autocorrelation is perfect, are set
                                    % to one (see corr(), lines 317-328). We
                                    % will replace these with zero,
                                    % although it is not technically
                                    % correct according to the above ref.
                                  
    else
        Rmat = RmatIn;
        Pmat = PmatIn;
        clear('RmatIn','PmatIn')
    end
    
    % Put 0's in place of 1's on Pmat diagonal (to preserve self-loops,
    % although pvals are not meaningful for these. Using linear indexing:
        Pmat(1:size(Pmat,1)+1:numel(Pmat)) = 0;
    
%% Later option: pass me an array of group columns with groups lining up with nodes or edges
    % call it additionalColumns
        % size = length(features.peaks.shifts), numCols
    % addnlColumnsNames = cell array of strings for each column
    % for each column:
        % varnames = 
                    for i = 1:size(additionalColumns,1)
                        eval(sprintf('%s = additionalColumns{i,2};',additionalColumns{i,1}));
                    end

    %% Generate the Network file
            % In case the columnLabels need to be a cell array of strings:
                %columnLabels = 1:size(features.intensities,2);
                %columnLabels = cellfun(@num2str,num2cell(features.peaks.shifts)','UniformOutput',0);
                %rowLabels = columnLabels;
            % For each element of edgesT > 0, list the columnLabel and the rowLabel as
                % Thresholding the edges
                        fprintf('\nApplying thresholding...\n')
                    switch whichCorr
                        case 'both'
                            % Get the indices for elements passing the
                                % cutoffs
                                passInd = ~or( abs(Rmat)<Rthresh,    Pmat>Pthresh   );
                            % Threshold Rmat
                                edgesR = Rmat;
                                edgesR(~passInd) = 0;
                                RmatT = edgesR;
                            % Threshold Pmat
                                edgesP = Pmat;
                                edgesP(~passInd) = 1;
                                PmatT = edgesP;
                            % Older, slower way:
                                %{
                                edgesR = Rmat;
                                edgesR(   or( abs(edgesR)<Rthresh,    Pmat>Pthresh   ) )= 0;
                                RmatT = edgesR;
                                edgesP = Pmat;
                                edgesP(   or( abs(edgesR)<Rthresh,    Pmat>Pthresh   ) )= 1;
                                PmatT = edgesP;
                            %}
                        case 'positive'
                            % Get the indices for elements passing the
                                % cutoffs
                                passInd = ~or( Rmat<Rthresh,    Pmat>Pthresh   );
                            % Threshold Rmat
                                edgesR = Rmat;
                                edgesR(~passInd) = 0;
                                RmatT = edgesR;
                            % Threshold Pmat
                                edgesP = Pmat;
                                edgesP(~passInd) = 1;
                                PmatT = edgesP;                            
                            % Older, slower way:
                                %{
                                edgesR = Rmat;
                                edgesR(   or( edgesR<Rthresh,    Pmat>Pthresh   ) )= 0;
                                RmatT = edgesR;
                                edgesP = Pmat;
                                edgesP(   or( edgesR<Rthresh,    Pmat>Pthresh   ) )= 1;
                                PmatT = edgesP;
                            %}
                        case 'negative'
                            % Get the indices for elements passing the
                                % cutoffs
                                passInd = ~or( (1-Rmat)<Rthresh,    Pmat>Pthresh   );
                            % Threshold Rmat
                                edgesR = Rmat;
                                edgesR(~passInd) = 0;
                                RmatT = edgesR;
                            % Threshold Pmat
                                edgesP = Pmat;
                                edgesP(~passInd) = 1;
                                PmatT = edgesP;                            
                            % Older, slower way:
                                %{
                                edgesR = Rmat;
                                edgesR(   or( (1-edgesR)<Rthresh,    Pmat>Pthresh   ) )= 0;
                                RmatT = edgesR;
                                edgesP = Pmat;
                                edgesP(   or( (1-edgesR)<Rthresh,    Pmat>Pthresh   ) )= 1;
                                PmatT = edgesP;
                            %}
                        otherwise
                            sprintf('\n\n\tIncorrect input received for variable "whichCorr". Selecting "both".\n\n')
                            % Get the indices for elements passing the
                                % cutoffs
                                passInd = ~or( abs(Rmat)<Rthresh,    Pmat>Pthresh   );
                            % Threshold Rmat
                                edgesR = Rmat;
                                edgesR(~passInd) = 0;
                                RmatT = edgesR;
                            % Threshold Pmat
                                edgesP = Pmat;
                                edgesP(~passInd) = 1;
                                PmatT = edgesP;
                            % Older, slower way:
                                %{
                                edgesR = Rmat;
                                edgesR(   or( abs(edgesR)<Rthresh,    Pmat>Pthresh   ) )= 0;
                                RmatT = edgesR;
                                edgesP = Pmat;
                                edgesP(   or( abs(edgesR)<Rthresh,    Pmat>Pthresh   ) )= 1;
                                PmatT = edgesP;
                            %}
                    end
                    
                fprintf('\nBuilding interactions table...\n')
                % Get the matrix row and column given the linear index of an element:
                    [row,col,Corr] = find(edgesR);
%                     Pval = edgesP(edgesP<1);
                    %CorrTemp = Rmat(passInd);
                    Pval = Pmat(passInd);
                    %[row,col] = ind2sub(size(Rmat),passInd);
                    Source = labels(col)';
                    Target = labels(row)';

                % Index Source and Target to facilitate sorting later on:
                    SourceInd = 1:size(Rmat,1);
                    TargetInd = 1:size(Rmat,2);      
                    SourceInd = SourceInd(col)';
                    TargetInd = TargetInd(row)';
                    EdgeName = [SourceInd TargetInd];
%                     addString = [];
%                     for i = 1:size(additionalColumns,1)
%                         addString = [addString,',',additionalColumns{i,1}];                        
%                     end
                     interactionsTable = table(Source,Target,Corr,Pval,EdgeName);
%                    eval(sprintf('interactionsTable = table(Source,Target,Corr,Pval,EdgeName%s);',addString));
                                   
                    % Remove the self edges (ignore direction)
%                     if strcmp(removeSelfEdges,'removeSelfEdges')
%                         interactionsTable(find(interactionsTable.Source == interactionsTable.Target),:) = [];
%                     end
% {
                    if strcmp(removeSelfEdges,'removeSelfEdges')                      
                        fprintf('\nRemoving Self Edges...\n')
                        interactionsTable(find(strcmp(interactionsTable.Source, interactionsTable.Target)),:) = [];
                    else 
                        fprintf('\nSkipping step to remove self edges...\n')
                    end
                    % Remove duplicate edges
                        % Sort edgename
                        fprintf('\nRemoving Duplicate Edges...\n')
                        interactionsTable.EdgeName = sort(interactionsTable.EdgeName,2);
                        [~,IA,~] = unique(interactionsTable.EdgeName,'rows');
                        interactionsTable.EdgeName = []; % clean this up
                        interactionsTable=interactionsTable(IA,:);      % remove duplicates
            %}            
toc

                % Convert table to cell array
                    %cellsTable = table2cell(interactionsTable);
                % Print the soft-thresholded network interactions to a .csv:
                    fprintf('\nWriting file...\n')
                    % For tab-delimited
                        writetable(interactionsTable,filename,'Delimiter','\t'); % tab-delimited
%                     % For csv
%                         writetable(interactionsTable,filename); % writes to .csv
%                             fprintf(['Generated ',num2str(size(interactionsTable,1)),' interactions\n\tin ',num2string(elapsedTime),' seconds.\n'])
%                             fprintf(['Network was written to "',filename '" in \n"',cd(),'"...\n\n']);
%                             clear('Source','Target','EdgeName','row','col','Corr','edges','edgesT','thresh')
%                     % For big tables
%                         cell2csv(filename,cellsTable)

                    %save('CNC_01_Preliminary_newFeatures_AlltoAll_to_NET.mat')
                    
                    end
