function arguments = parse_args(defaults, varargin)

% Example: 
% defaults = struct('minPeakHeight',5000,'c12AR',[.02 .06],'verbose',0,'normalizeMethod','logMR');
% args = parse_args(defaults, {'verbose', 1, 'norm', 'total'})

% Parse input options
optionNames = fieldnames(defaults);
if ~isempty(varargin)
    for pair = reshape(varargin{:},2,[]) % pair is {propName;propValue}
        inpName = lower(pair{1}); % make case insensitive
        m=strcmpi(validatestring(inpName,optionNames),optionNames); %allows partial unambiguous matches
        if any(m)
            defaults.(optionNames{m})=pair{2};
        else
            error('%s is not a recognized parameter name',inpName)
        end
    end
end
arguments = defaults;