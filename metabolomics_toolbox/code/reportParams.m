function p = reportParams(varargin)
%% reportParams

% This function, when called at the beginning of a function, deals the
% input parameters to an output structure, p. This allows us to record
% parameters effectively and dynamically, thereby facilitating
% reproducibility in workflows. This function is intended to be used before
% anything happens in a function so that inputs can be accurately captured,
% but may be placed after parsing varargins and assigning defaults. If
% pasted before varargins are parsed, only the param 'varargin' will be 
% recorded (as a cell array). 
% 
% Note: In order to avoid making lots of copies of large matrices, a
% sizeLimit option is provided. Default is 10000 elements per input param.
% If the numel() of the input param exceeds the limit, a string is returned
% for that param within the p structure indicating that the data were too
% large. This is particularly useful when running a function within a loop.
% 
% Note: This function may cause a reduction in speed. This is difficult to
% test as it depends on the number of vars and their sizes. Likewise, be
% mindful of memory useage. 
%
% Optional Params (name-value pair): 
%           'sizeLimit'     numeric (double, int, etc.) specifying the 
%                           maximum number of elements allowed to be stored 
%                           for each param. ** If cell array of doubles,
%                           only the number of cells is counted** 
%           'exclude'       cell array of strings specifying the names of 
%                           parameters to exclude in the output (e.g. if
%                           they are too big, or are meant to be inputs
%                           rather than parameters). 
% Usage:
%
%     p = reportParams();                                    % using default sizeLimit
%     p = reportParams('sizeLimit',3);                       % sizeLimit 3 elements
%     p = reportParams('exclude',{'X','ppm'});               % exclusion
%     p = reportParams('exclude',{'X','ppm'},'sizeLimit',3); % combined is fine
%                                                             
%
% MTJ 2020


%% Parse Name-value pairs (if present)

    if ~isempty(varargin)
        for a = 1:2:length(varargin)
            switch varargin{a}
                
                case 'sizeLimit'
                    sizeLimit = varargin{a+1};
                    
                case 'exclude'
                    excludeParams = varargin{a+1};
                    
                % other optional params can be added later as case
                % statements
            end
        end
    end
        
   if ~exist('sizeLimit','var')
       sizeLimit = 10000;
   end
    
%%   
    params = evalin('caller','who');
    if exist('excludeParams','var')
        params = params(~ismember(params,excludeParams));
    end
    p = struct(); % default is empty struct
    
    if numel(params)>0
    %% Run this to copy all params and corresponding values to an output struct 'p'

        for i = 1:length(params)
           p.(params{i}) = evalin('caller',params{i});
        end

    %% Run this to avoid making another copy of ~large data
        for i = 1:length(fieldnames(p))
            if numel(p.(params{i})) > sizeLimit
                p.(params{i}) = ['> ',num2str(sizeLimit),' elements. Data not stored.'];            
            end
        end
    
    end
    
end