%% Wormsorter workflow
% Francesca V. Ponce 
% June 2016
% ->PROJECT ID HERE<-
%
% *** After adding your project ID, save a copy of this file into your data
% analysis folder****
%
%***After running script, save workspace variables in your data analysis
%folder****

%%
% selects the folder the .txt sorter files are saved in
% # build a list of file names with absolute path
fPath = uigetdir('.', 'Select directory');
if fPath==0, error('no folder selected'), 
end
fNames = dir( fullfile(fPath,'*.txt') );
filenames = strcat(fPath, filesep, {fNames.name});
filenames = filenames';
shortfilenames = {fNames.name};
shortfilenames = shortfilenames';

%% 
% 'parameter_list.mat' is a .mat file with the list of import functions of 
% parameters to be imported

load('parameter_list.mat');
parameters = cell(length(filenames),1);
for n = 1:length(filenames);
    parameters{n,1} = import_parameters(filenames{n,1},parameter_list);
end

%%
% getting instrument settings of each file and storing in a matrix
%[sorter_settings] = get_settings(parameters);
[sorter_settings, ext_detector_power] = get_settings(parameters);
% erasing blank cells
[parameters,numericparameters,non_numericparameters, numericparameters_list] = check_parameters(parameters,parameter_list);

%%
% Getting control particle data and generate plots of control particle data 
%  **If you don't want to see the plots, change the number "1" at the end
%  of the argument to a "0"
[controlparticles, controlparticles_parameters_list, mean_controlparticles, stds_controlparticles, intersectt] = findcontrolparticles...
    (numericparameters, numericparameters_list, sorter_settings, 1);
%%

