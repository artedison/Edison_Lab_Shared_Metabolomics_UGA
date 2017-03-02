function [parameters,numericparameters,non_numericparameters,numericparameters_list] = check_parameters(parameters,parameter_list)

% checks which parameters are numeric and which ones are not.

numericparameters =parameters;
non_numericparameters = parameters;
numericparameters_list = parameter_list;

for j = 1:length(parameters);
    parametersna{j,1} = cellfun(@isnumeric,(parameters{j,1}));
    parametersna_zeroes{j,1} = find(~parametersna{j,1});
    parametersna_ones{j,1} = find(parametersna{j,1});
    numericparameters{j,1}(parametersna_zeroes{j,1})=[];
    non_numericparameters{j,1}(parametersna_ones{j,1})=[];
end

% erases empty cells at the bottom of the array of numericparametes
numericparameters_list(parametersna_zeroes{j,1}) = [];
fh = @(x) all(isnan(x(:)));
numberofnumericparameters = numericparameters{1,1};
for i = 1:1:length(numericparameters);
    for k = 1:length(numberofnumericparameters);
        numericparameters{i,1}{k,1} = num2cell(numericparameters{i,1}{k,1});
        numericparameters{i,1}{k,1}(cellfun(fh, numericparameters{i,1}{k,1})) = [];
    end
end