function [relativePath] = fullPath2RelativePath(startPath_full,endPath_full,varargin)
%% fullPath2RelativePath
% 
%   Produce the relative path string necessary to navigate from one 
%   directory to another. Optional flag argument to convert to a
%   shell-friendly string by adding escape characters when necessary. This
%   is particularly useful when passing paths as arguments to shell scripts
%   from within Matlab, and eliminates common human error. Furthermore, it
%   allows relative paths to be used within workflows where directory
%   structures may change over time, since relative paths are re-calculated
%   during each run. 
%
%   MTJ DEC2020
%
%   NOTE: both dirs must actually have a common root
%   Potential Bugs:
    % What happens if the locations are in the same dir? 
    % What happens if the one dir is within another?
    
%%
    % StrSplit both full paths based on / character
        spf_split = strsplit(startPath_full,'/');
            spf_split = spf_split(~cellfun(@isempty,spf_split));
        
        epf_split = strsplit(endPath_full,'/');
            epf_split = epf_split(~cellfun(@isempty,epf_split));
            
    % Starting from root, find the fork point (deepest common directory)
        commonPath = false(1,min([length(spf_split),length(epf_split)]));
        
        for i = 1:length(commonPath)
            commonPath(i) = strcmp(spf_split{i},epf_split{i});
        end
        
        % Unbroken true
             forkInd = find(commonPath == 0,1,'first');
        
    % Using the index of the fork dir, calculate the number of pops to get
    % there from starting location. Then, use the directory names on the
    % end location side to navigate from that point, and build the relative 
    % path string:
        partialPath_end = join(epf_split(forkInd:end),'/');
        partialPath_start = repmat('../',1,length(spf_split(forkInd:end)));
        relativePath = [partialPath_start,...
                        partialPath_end{:}];
                    
    % Optional shell scripting conversion
        if contains(varargin,'useEscapeCharacters')
            relativePath = regexprep(relativePath,' ','\\ ');
            relativePath = regexprep(relativePath,'(','\\(');
            relativePath = regexprep(relativePath,')','\\)');
        end
                    
end