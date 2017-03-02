function [data] = deletecontrolparticles(data)
load('intersecttest.mat');
for z=1:length(intersectt);
    inter = intersectt{z,1}';
    data(:,inter) = []; 
end