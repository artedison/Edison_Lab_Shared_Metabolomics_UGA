function parameters = import_parameters(filename,parameter_list)

% imports the parameters listed in parameters_list
% each column in 'parameters' contains values of each parameter

parameters = cell(length(parameter_list),1);
for i = 1:length(parameter_list);
    sh = str2func(parameter_list{i,1});
    parameters{i,1} = sh(filename);
end