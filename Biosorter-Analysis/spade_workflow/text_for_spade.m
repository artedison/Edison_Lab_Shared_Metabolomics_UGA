%% selects the folder the .txt files are saved in
%# build a list of file names with absolute path
fPath = uigetdir('.', 'Select directory');
if fPath==0, error('no folder selected'), 
end
fNames = dir( fullfile(fPath,'*.txt') );
filenames = strcat(fPath, filesep, {fNames.name});
shortfilenames = {fNames.name};
load('intersectt.mat');
intersectt = intersectt';
%%
for i= 1:length(filenames);
    filename = filenames{i};
    [filename, data, marker_names,channel_names] = textfilestomat(filename);
    inter = intersectt{1,i};
    data(:,inter) = [];
    filename1 = sprintf('sorterfile%d', i);
    save(filename1, 'filename','data','marker_names','channel_names');
    writefcs(filename, data, marker_names,channel_names);
end
%%
SPADE
%%