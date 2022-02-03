function [spectra,matchedMetadata,combinedStruct] = mapMetadataToSpectra(spectra,mdTable,mdColumn)

%% Map a metadata table to the spectra struct by experiment number
%   NOTE: The metadata table mdTable may contain more rows than the
%           spectra struct. Only the matching rows will be used. Uses the
%           'Title' field in spectra to match a field in 'mdTable'. mdTable
%           fields cannot overlap with fields in 'spectra'. 
%
%   Inputs:
%         spectra = spectra;        % spectra structure (e.g. the result of
%                                   % loadallft or Load1D). Field 'Title'
%                                     required; must match the provided
%                                     metadata field contents
%         mdTable = metadataTable;  % metadataTable (must include runIDs)
%         mdfield = 'projectID';    % In the case of several runIDs, the
%                                   % name of the one to be used for
%                                   % mapping (variable name in mdTable)
%   Outputs:
%         combinedStruct            a struct array containing with the metadata for matching 
%                                   rows concatenated to the spectra struct array.
%   Usage:
%         combinedStruct = [combinedStruct] = mapMetadataToSpectra(spectra,mdTable,'projectID');
%   NOTE: 
% MTJ 2020

%% 
            % Add titles as doubles into the spectra struct.
                tAd = num2cell(cellfun(@str2double, {spectra.Title}')); % convert cells with strings to cells with doubles
                [spectra(:).titlesAsDoubles] = tAd{:};                  % make these a new field in spectra struct
                clear('tAd')
            % Find the matches between projectID and spectra titles    
                [lia,locb] = ismember([spectra.titlesAsDoubles],mdTable.(mdColumn)); % *** EXPERIMENT-SPECIFIC ***
            
                % Report on whether or not all were matched; indicate those
                % which were not.
                    if ~all(lia)
                        fprintf(['\n\tThe following spectra did not have a match:\n\t\t',...
                                num2str([spectra(~lia).titlesAsDoubles]),'\n'])
                    else
                        fprintf('\n\tAll spectra were matched\n\n')
                    end
                    
            % Sort them both according to experiment ID (just to be sure)
                combinedStruct = table2struct([mdTable(locb(lia),:),struct2table(spectra(lia))]);
                matchedMetadata = mdTable(locb(lia),:);
                
end