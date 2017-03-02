%% Wormsorter workflow
% Francesca V. Ponce 
% June 2016
% ->PROJECT ID HERE<-
%
% After adding your project ID, save a copy of this file into your data
% analysis folder then you can modify the workflow as needed.
% Analysis pipelines for the data will change from project to project and 
% the pipeline described here may not be the best fit for your data,
% but it can be adapted to your data.
%%
% selects the folder the .txt sorter files are saved in
% build a list of file names with absolute path
% these are the files that will be analyzed
fPath = uigetdir('.', 'Select directory');
if fPath==0, error('no folder selected'), 
end
fNames = dir( fullfile(fPath,'*.txt') );
filenames = strcat(fPath, filesep, {fNames.name}');
shortfilenames = {fNames.name}';
%%
% getting instrument settings of each file and storing in a matrix
[sorter_settings] = get_sortersettings(filenames);
%% 
% 'parameter_list.mat' is a .mat file with the list of import functions of 
% parameters to be imported
load('parameter_list.mat');
parameters = cell(length(filenames),1);
for n = 1:length(filenames);
    parameters{n,1} = import_parameters(filenames{n,1},parameter_list);
end
%%
% erasing blank cells and getting numeric and non numeric parameters
[parameters,numericparameters,non_numericparameters, numericparameters_list] = check_parameters(parameters,parameter_list);

%%
% Getting control particle data and generate plots of control particle data 
% If you don't want to see the control particle plots, delete the number "1" at the end
% of the input argument
[controlparticles, controlparticles_parameters_list, mean_controlparticles, stds_controlparticles, intersectt] = findcontrolparticles...
    (numericparameters, numericparameters_list, sorter_settings, 1);

%% plots raw data points according to plot_rawdata function
raw_data_plot = input('Plot raw data (Y/N)?', 's');
if  strcmp(raw_data_plot,'Y');
    plot_rawdata(parameters, filenames)
end
%%
