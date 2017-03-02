function plot_controlparticles(controlparticles, sorter_settings, intersectt, mean_controlparticles, stds_controlparticles)

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
