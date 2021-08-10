function [repInds,uniqueRows,mdSubset,included] = indexReps(mdTable,variableNames,varargin)

%% indexReps

% Groups replicates, based on desired fields. 
%   
%       * only exact matches are supported
%
% Inputs: 
%         mdTable         metadata table containing variableName
%         variableNames   names of the variables (column) in mdTable by
%                         which replicates are defined. Must be of the
%                         form: 
%                               {'VariableName1','VariableName2'}
%         varargin options:        
%             'exclude'   specific names of values in mdTable.variableName 
%                         to be excluded Must be in name-value pair form as
%                         follows:
%                               'exclude',{'VariableName','StringToExclude'}
%                               'exclude',{'VariableName',numberToExclude}
%                               'exclude',{'VariableName',logicalToExclude}
%
%                         NOTE: excluded option can be layered (see usage
%                         below)
%                         NOTE: function may break if exclude option
%                         eliminates a column specified by variableName
%             'only'      the inverse of 'exclude'. Only includes those
%                         rows for which the specified variable matches the
%                         value provided. This can also be layered, but do
%                         so with caution. With multiple inputs,
%                         intracategory onlys are interpreted as OR, while
%                         intercategory onlys are interpreted as AND.
% 
% Outputs:
%         repInds         vector of group indices for each row in mdTable.
%                         Any excluded (not grouped) rows are represented
%                         as zeros.
%         uniqueRows      mdTable subset with one row for each unique group
%                         including repIndex for each group as used in
%                         repInds.
%
% 
% Usage:
%
%     [repInds,uniqueRows,mdSubset,included] = indexReps(md,{'Genotype','Aeration'},...
%                                               'exclude',{'Genotype','Pooled'},...
%                                               'exclude',{'CollectionState','D'},...
%                                               'exclude',{'CollectionState','NA'});
% 
%     [repInds,uniqueRows,mdSubset,included] = indexReps(md,{'Genotype'},...
%                                               'only',{'Genotype','bd'},...
%                                               'exclude',{'Genotype','Media Blank'},...
%                                               'exclude',{'Aeration','Sealed'},...
%                                               'exclude',{'CollectionState','D'},...
%                                               'exclude',{'CollectionState','NA'});
% 
%     
% MTJ 2020



%% Parse varargin options

        if ~isempty(varargin)        
            for i = 1:2:length(varargin)

                switch varargin{i}
                    case 'only'                       
                        only{i,:} = varargin{i+1};
                        
                        %pull = castToCells(varargin{i+1});
                        %pull = varargin{i+1};
                        only(cellfun(@isempty,only)) = []; % clean up empty cells
                        
                    case 'exclude'
                        exclude(i) = varargin(i+1);
                        %exclude = castToCells(varargin{i+1});
                        %exclude = varargin{i+1};
                        exclude(cellfun(@isempty,exclude)) = []; % clean up empty cells
                end
            end
            
        end             


%% Exclude samples if necessary (decide which rows to index)

    % After running this, 'include' will always be set
        
        if and(exist('exclude','var'),~isempty('exclude')) % if excluding
            
            % Loop through exclude args and note the rows that are exluded
                excluded = zeros(height(mdTable),length(exclude));
            
                for i = 1:length(exclude)
                    [~,excluded(:,i)] = ismember(mdTable.(exclude{i}{1}),exclude{i}{2});
                end
            % Convert to included rows
                notExcluded = ~any(excluded,2); % want all rows that weren't excluded by ANY 'exclude' condition
                
                
        else    % if nothing to exclude, set to all rows:
                exclude = {};
                notExcluded = ones(height(mdTable),1);
        end

%% Exclude samples if 'only' is used (decide which rows to index)

    % After running this, 'onlys' will always be set
        
        if and(exist('only','var'),~isempty('only')) % if excluding
            
            % Divide only inputs into fields they're applied to 
                onlyIns = [only{:}];
                [onlyfields,~,ofInds] = unique(onlyIns(1:2:end));
            
            % We'll require AND across categories, but OR within them
                metaInclude = zeros(height(mdTable),length(onlyfields));
            
            for c = 1:length(onlyfields) % for each category of only inputs
                
                catInds = find(ofInds == c); % inds for the args pertaining to current category 
                
                % Loop through args and note the rows that are exluded
                    include = zeros(height(mdTable),numel(catInds));

                    for i = 1:length(catInds)
                        [~,include(:,i)] = ismember(mdTable.(only{catInds(i)}{1}),only{catInds(i)}{2});
                    end

                % Convert to included rows       
                    metaInclude(:,c) = any(include,2); % must satisfy at least one 'only' condition (need OR for intracategory onlys; AND for intercategory onlys)
                
            end
            
            % Now apply AND across only categories:
            
                onlys = all(metaInclude,2);
            
        else    % if nothing to exclude, set to all rows:
                only = {};
                onlys = ones(height(mdTable),1); 
        end
        
%% Combine 'exclude' and 'only' 

    included = and(notExcluded,onlys);
        
%% Index the rows 
    
    % Find the inds for the columns we want (indicated in variableNames)
        allVars = mdTable.Properties.VariableNames;
        [~,varInds] = ismember(variableNames,allVars);
        
        % Also get the inds for those vars in 'onlys' and 'exclude' which 
        % do not overlap with variableNames
            onlyDummy = [only{:}];
            excludeDummy = [exclude{:}];
            otherVars = [reshape(onlyDummy,2,[]),reshape(excludeDummy,2,[])];
            if ~isempty(otherVars)
                [~,otherVarInds] = ismember(setdiff(unique(otherVars(1,:)),variableNames),allVars);
            else
                otherVarInds = [];
            end
    % ID the reps for the included rows, set rest to zeros
        [uniqueRows,~,a] = unique(mdTable( included ,varInds),'rows');
        %[uniqueRows,~,a] = unique(mdTable( included ,[varInds(:);otherVarInds(:)]),'rows');
        
        uniqueRows.RepIndex = (1:height(uniqueRows))';
        
        repInds = double(included);
        repInds(included) = a; 
        
    % Add in other fields that were mentioned
        
        mdSubset = mdTable(included ,[varInds(:);otherVarInds(:)]);
        mdSubset.RepIndex = repInds(repInds>0);
        mdSubset.OriginalRows = find(repInds>0);
    
end
