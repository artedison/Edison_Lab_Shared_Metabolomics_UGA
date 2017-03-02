function [filename, data, marker_names,channel_names] = textfilestomat(filenames)

fh = @(x) all(isnan(x(:)));
    
filename = filenames;
    
TOF = import_tof(filename);
Extinction = import_ext(filename);
Green = import_green(filename);
Red = import_red(filename);
Yellow = import_yellow(filename);
PHExtinction = import_PHExtinction(filename);
PWExtinction = import_PWExtinction(filename);
PHGreen = import_PHGreen(filename);
PWGreen = import_PWGreen(filename);
PHRed = import_PHRed(filename); 
PWRed = import_PWRed(filename);
PHYellow = import_PHYellow(filename);
PWYellow = import_PWYellow(filename);
    
    

data = {TOF';Extinction';Green';Red';Yellow';PHExtinction';PWExtinction';PHGreen';PWGreen';PHRed';PWRed';PHYellow';PWYellow'};

% if you want to take log of data, uncomment line 28

    size_data = size(data);
    for i= 1:size_data(1,1);
    data{i,1} = num2cell(data{i,1});
    data{i,1}(cellfun(fh, data{i,1})) = [];
    data{i,1} = cell2mat(data{i,1});
    %data{i,1} = log10(data{i,1});    
    end
    data = cell2mat(data);
    
    marker_names = {'TOF'; 'Extinction';'Green';'Red';'Yellow';'PHExtinction';'PWExtinction';'PHGreen';'PWGreen';'PHRed';'PWRed';'PHYellow';'PWYellow'};
    channel_names = {'TOF'; 'Extinction';'Green';'Red';'Yellow';'PHExtinction';'PWExtinction';'PHGreen';'PWGreen';'PHRed';'PWRed';'PHYellow';'PWYellow'};
 
 