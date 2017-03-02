%% Binning workflow
% Francesca V. Ponce
% July 2016
% PROJECT ID HERE
%% Be sure that workspace variables from wormstorter_workflow are loaded before running this script

%%
% ext = cell array raw extinction values from data files with gaps from control
% particles
% Ext = mat variable with empty spaces erased
% nExt = normalized extinction values
%nExt
%% 
%create new matrix of events without control particles
numericparameters_nocp = numericparameters;
for l = 1:length(numericparameters_nocp);
    for k = 1:length(numericparameters_nocp{l,1});
    numericparameters_nocp{l,1}{k,1}([intersectt{l,1}],:) = []; 
    end
end
%%
% taking extinction values from  numericparameters_nocp and putting them in
% new matrix (ext)
for l = 1:length(numericparameters_nocp); 
   ext{l,1} = numericparameters_nocp{l,1}{7,1}; 
end
%%
%taking TOF values from numericparameters_nocp into new matrix (tof)
for l = 1:length(numericparameters_nocp); 
   tof{l,1} = numericparameters_nocp{l,1}{6,1}; 
end

%%
% manipulating extinction value matrix for processing w prexisting code
% blocks
ext = ext';

for j = 1:length(ext);
Ext{1,j} = cell2mat(ext{1,j});
end
%%
%take mean extinction values from mean_controlparticles into new double
%with correct dimensions
mean_controlparticles_ext = mean_controlparticles(:,2);
mean_controlparticles_ext = mean_controlparticles_ext';

%% Normalization of extinction values by dividing all extinctions values in each file by the mean control particle extinction from that same file
for h = 1:length(mean_controlparticles_ext);
    Ext{1,h} = num2cell(Ext{1,h});
    nExt{1,h} = cellfun(@(v) v./mean_controlparticles_ext{1,h}, Ext{1,h},'UniformOutput', false);
end

%% take log of normalized extinction values
for l = 1:length(nExt);
   nExt{1,l} = cell2mat(nExt{1,l});
   lnExt{1,l} = log(nExt{1,l});    
end

%% invoke worm binning script to generate table with number of events in each defined bin

make_worm_bins

%% plotting histograms
[rows columns] = size(m);
for k = 1:rows;
    subplot(1,rows,k)
    bar(m(k,:));
end