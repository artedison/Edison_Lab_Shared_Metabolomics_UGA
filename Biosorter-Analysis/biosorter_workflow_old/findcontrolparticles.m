function [controlparticles, controlparticles_parameters_list, mean_controlparticles, stds_controlparticles, intersectt] = findcontrolparticles(numericparameters, numericparameters_list, sorter_settings, makeplot)
% finds control particles
% input: matrix of numeric parameters of worm sorter data
% output: 'control particles' matrix of control particles and
% 'info_controlparticles' with the mean, standard deviation of the control
% particles in each file

% find the parameters that will be used to detect control particles

index_TOF = strfind(numericparameters_list, 'import_TOF');
index_Extinction = strfind(numericparameters_list, 'import_Extinction');
index_Green= strfind(numericparameters_list, 'import_Green');
index_Yellow = strfind(numericparameters_list, 'import_Yellow');
index_Red = strfind(numericparameters_list, 'import_Red');
index_PHExtinction = strfind(numericparameters_list,'import_PHExtinction');
index_PWExtinction = strfind(numericparameters_list,'import_PWExtinction');
index_PHGreen = strfind(numericparameters_list,'import_PHGreen');
index_PWGreen = strfind(numericparameters_list,'import_PWGreen');
index_PHYellow = strfind(numericparameters_list,'import_PHYellow');
index_PWYellow = strfind(numericparameters_list,'import_PWYellow');
index_PHRed = strfind(numericparameters_list,'import_PHRed');
index_PWRed = strfind(numericparameters_list,'import_PWRed');

index_TOF = find(~cellfun(@isempty,index_TOF));
index_Extinction = find(~cellfun(@isempty,index_Extinction));
index_Green = find(~cellfun(@isempty,index_Green));
index_Yellow = find(~cellfun(@isempty,index_Yellow));
index_Red = find(~cellfun(@isempty,index_Red));
index_PHExtinction = find(~cellfun(@isempty,index_PHExtinction));
index_PWExtinction = find(~cellfun(@isempty,index_PWExtinction));
index_PHGreen = find(~cellfun(@isempty,index_PHGreen));
index_PWGreen = find(~cellfun(@isempty,index_PWGreen));
index_PHYellow = find(~cellfun(@isempty,index_PHYellow));
index_PWYellow = find(~cellfun(@isempty,index_PWYellow));
index_PHRed = find(~cellfun(@isempty,index_PHRed));
index_PWRed = find(~cellfun(@isempty,index_PWRed));

index_controlparticles_parameters = [index_TOF; index_Extinction; index_Green; ...
    index_Yellow; index_Red; index_PHExtinction; index_PWExtinction; index_PHGreen; ...
    index_PWGreen; index_PHYellow; index_PWYellow; index_PHRed; index_PWRed];
 
% list of parameters that are needed to identify control particles

controlparticles_parameters_list = {'TOF'; 'Extinction'; 'Green'; ...
    'Yellow'; 'Red'; 'PHExtinction'; 'PWExtinction'; 'PHGreen'; ...
    'PWGreen'; 'PHYellow'; 'PWYellow'; 'PHRed'; 'PWRed'};

% make a matrix with the parameters that will be used to detect control particles

for i = 1:length(numericparameters);
    for k = 1:length(index_controlparticles_parameters);
        controlparticles{i,1}{k,1} = numericparameters{i,1}{(index_controlparticles_parameters(k)),1};
        controlparticles{i,1}{k,1} = cell2mat(controlparticles{i,1}{k,1});
    end
end

% thresholds to detect control particles

for i = 1:length(controlparticles);
    tindx{1,i} = find((controlparticles{i,1}{1,1}>50)); %  TOF
    eindx1{1,i} = find((controlparticles{i,1}{2,1}>10)); % Extinction
    eindx{1,i} = find((controlparticles{i,1}{2,1}<60));
    gindx{1,i} = find((controlparticles{i,1}{3,1}>10)); % Green
    yindx{1,i} = find((controlparticles{i,1}{4,1}>10)); % Yellow
    rindx{1,i} = find((controlparticles{i,1}{5,1}>10)); % Red
    ephindx{1,i} = find((controlparticles{i,1}{6,1}>3000)); % PHExtinction
    gphindx{1,i} = find((controlparticles{i,1}{8,1}>3000)); % PHGreen
    yphindx{1,i} = find((controlparticles{i,1}{10,1}>3000)); % PHYellow
    rphindx{1,i} = find((controlparticles{i,1}{12,1}>3000)); % PHRed
end

% making all vectors the same size to intersect
% length(controlparticles) would be the numer of files 

for k = 1:length(controlparticles);
    
    if length(eindx1{1,k})<length(eindx{1,k});
    eindx1{1,k}(numel(eindx{1,k})) = 0;
    elseif length(eindx{1,k})<length(eindx1{1,k});
    eindx{1,k}(numel(eindx1{1,k})) = 0;
    end
    
    if length(eindx{1,k})<length(gindx{1,k});
    eindx{1,k}(numel(gindx{1,k})) = 0;
    elseif length(gindx{1,k})<length(eindx{1,k});
    gindx{1,k}(numel(eindx{1,k})) = 0;
    end
    
    if length(eindx{1,k})<length(yindx{1,k});
    eindx{1,k}(numel(yindx{1,k})) = 0;
    elseif length(yindx{1,k})<length(eindx{1,k});
    yindx{1,k}(numel(eindx{1,k})) = 0;
    end
    
    if length(eindx{1,k})<length(rindx{1,k});
    eindx{1,k}(numel(rindx{1,k})) = 0;
    elseif length(rindx{1,k})<length(eindx{1,k});
    rindx{1,k}(numel(eindx{1,k})) = 0;
    end
    
    if length(eindx{1,k})<length(gphindx{1,k});
    eindx{1,k}(numel(gphindx{1,k})) = 0;
    elseif length(gphindx{1,k})<length(eindx{1,k});
    gphindx{1,k}(numel(eindx{1,k})) = 0;
    end
    
end

% intersect to find the common indices

for i = 1:length(controlparticles);
intersect1{1,i} = intersect(eindx{1,i},gindx{1,i});
intersect11{1,i} = intersect(eindx1{1,i},intersect1{1,i});
intersect2{1,i} = intersect(intersect11{1,i},yindx{1,i});
intersect3{1,i} = intersect(intersect2{1,i},rindx{1,i});
intersect4{1,i} = intersect(intersect3{1,i},gphindx{1,i});
intersectt{i,1} = intersect(intersect3{1,i},intersect4{1,i});
end

% finding controlparticles

for i = 1:length(controlparticles);
    for k = 1:length(index_controlparticles_parameters);
        controlparticles{i,1}{k,1} = controlparticles{i,1}{k,1}(intersectt{i,1});
    end
end


for i = 1:length(controlparticles);
    for k = 1:length(index_controlparticles_parameters);
        mean_controlparticles(i,k) = mean(controlparticles{i,1}{k,1});
        stds_controlparticles(i,k) = std(controlparticles{i,1}{k,1});
    end
end
mean_controlparticles = num2cell(mean_controlparticles);
stds_controlparticles = num2cell(stds_controlparticles);

% controlparticles_summary_headers = {'TOF', 'Extinction', 'Green', ...
%     'Yellow', 'Red', 'PHExtinction', 'PWExtinction', 'PHGreen', ...
%     'PWGreen', 'PHYellow', 'PWYellow', 'PHRed', 'PWRed'};

mean_controlparticles = [controlparticles_parameters_list';mean_controlparticles];
stds_controlparticles = [controlparticles_parameters_list';stds_controlparticles];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ploting

if exist('makeplot')

sorter_settings = sorter_settings(2:end,:);
ext_detector_power = sorter_settings(:,1);
mean_controlparticles = mean_controlparticles(2:end,:);
stds_controlparticles = stds_controlparticles(2:end,:);

figure()

subplot(3,2,1)
for k = 1:length(controlparticles);
plot(controlparticles{k,1}{1,1},controlparticles{k,1}{2,1},'bo', 'color', rand(1,3));
hold on
xlabel('TOF')
ylabel('Ext')
end

subplot(3,2,2)
mean_TOF_cp = cell2mat(mean_controlparticles(:,1));
std_TOF_cp = cell2mat(stds_controlparticles(:,1));
bar(mean_TOF_cp, 'black');
hold on
errorbar(mean_TOF_cp,std_TOF_cp,'bo','Color', 'black')
hold on
xlabel('samples')
ylabel('mean TOF')

subplot(3,2,3)
mean_Extinction_cp = cell2mat(mean_controlparticles(:,2));
std_Extinction_cp = cell2mat(stds_controlparticles(:,2));
bar(mean_Extinction_cp, 'blue');
hold on
errorbar(mean_Extinction_cp,std_Extinction_cp,'bo','Color', 'blue')
hold on
xlabel('samples')
ylabel('mean Extinction')

subplot(3,2,4)
mean_Green_cp = cell2mat(mean_controlparticles(:,3));
std_Green_cp = cell2mat(stds_controlparticles(:,3));
bar(mean_Green_cp, 'green');
hold on
errorbar(mean_Green_cp,std_Green_cp,'bo','Color', 'green')
hold on
xlabel('samples')
ylabel('mean Green')

subplot(3,2,5)
mean_Yellow_cp = cell2mat(mean_controlparticles(:,4));
std_Yellow_cp = cell2mat(stds_controlparticles(:,4));
bar(mean_Yellow_cp, 'yellow');
hold on
errorbar(mean_Yellow_cp,std_Yellow_cp,'bo','Color', 'yellow')
hold on
xlabel('samples')
ylabel('mean Yellow')

subplot(3,2,6)
mean_Red_cp = cell2mat(mean_controlparticles(:,5));
std_Red_cp = cell2mat(stds_controlparticles(:,5));
bar(mean_Red_cp, 'red');
hold on
errorbar(mean_Red_cp,std_Red_cp,'bo','Color', 'red')
hold on
xlabel('samples')
ylabel('mean Red')


% controlparticles figures
figure()

subplot(3,2,1)
mean_PHExtinction_cp = cell2mat(mean_controlparticles(:,6));
std_PHExtinction_cp = cell2mat(stds_controlparticles(:,6));
bar(mean_PHExtinction_cp, 'blue');
hold on
errorbar(mean_PHExtinction_cp,std_PHExtinction_cp,'bo','Color', 'blue')
hold on
xlabel('samples')
ylabel('mean PHExtinction')

subplot(3,2,2)
mean_PHGreen_cp = cell2mat(mean_controlparticles(:,7));
std_PHGreen_cp = cell2mat(stds_controlparticles(:,7));
bar(mean_PHGreen_cp, 'green');
hold on
errorbar(mean_PHGreen_cp,std_PHGreen_cp,'bo','Color', 'green')
hold on
xlabel('samples')
ylabel('mean PHGreen')

subplot(3,2,3)
mean_PHYellow_cp = cell2mat(mean_controlparticles(:,8));
std_PHYellow_cp = cell2mat(stds_controlparticles(:,8));
bar(mean_PHYellow_cp, 'yellow');
hold on
errorbar(mean_PHYellow_cp,std_PHYellow_cp,'bo','Color', 'yellow')
hold on
xlabel('samples')
ylabel('mean PHYellow')

subplot(3,2,4)
mean_PHRed_cp = cell2mat(mean_controlparticles(:,9));
std_PHRed_cp = cell2mat(stds_controlparticles(:,9));
bar(mean_PHRed_cp, 'red');
hold on
errorbar(mean_PHRed_cp,std_PHRed_cp,'bo','Color', 'red')
hold on
xlabel('samples')
ylabel('mean PHRed')

subplot(3,2,5)
extdetectorpower = cell2mat(ext_detector_power);
plot(extdetectorpower,mean_Extinction_cp,'bo','Color', 'black');
xlabel('extinction detector power')
ylabel('mean Extinction')

subplot(3,2,6)
for k = 1:length(controlparticles);
plot(k,length(intersectt{k,1}),'bo','Color', 'black');
hold on
xlabel('samples')
ylabel('mean Extinction')
end

else
end