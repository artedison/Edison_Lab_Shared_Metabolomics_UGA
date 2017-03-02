function [parameters,numericparameters,non_numericparameters] = verify_parameters(parameters)

% checks which paramters are numeric and which ones are not.

numericparameters =parameters;
non_numericparameters = parameters;

for j = 1:length(parameters);
    parametersna{j,1} = cellfun(@isnumeric,(parameters{j,1}));
    parametersna_zeroes{j,1} = find(~parametersna{j,1});
    parametersna_ones{j,1} = find(parametersna{j,1});
    numericparameters{j,1}(parametersna_zeroes{j,1})=[];
    non_numericparameters{j,1}(numericparameters{j,1})=[];
end

