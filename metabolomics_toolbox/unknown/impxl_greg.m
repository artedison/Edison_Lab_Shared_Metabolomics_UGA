%% Import data from Peak pick worksheet

path{1}='N2_INAD1';
path{2}='N2_INAD2';
path{3}='N2_INAD3';
path{4}='N2_INAD4';
path{5}='N2_HS_INAD1';
path{6}='N2_HS_INAD2';
path{7}='N2_HS_INAD3';
path{8}='N2_HS_INAD4';
%%
for i=1:8
[~, ~, raw] = xlsread('C:\Users\Camera\Dropbox\Greg-Chaevien\Chae_problem1\lightintense.xlsx','N2_HS_INAD4');
spec = raw(2:end,1:3);
% make output 
data(:,:) = reshape([spec{:}],size(spec));

% array to column variable names
ppm1 = data(:,2);
DQ = data(:,3);

% Make Spectra structure

spectra(i).ppm1=ppm1;
spectra(i).ppm2=DQ;
Intensity=zeros(400,400);
Intensity(sub2ind(size(Intensity),round(ppm1), round(DQ)))=data(:,1);
spectra(i).real=Intensity;
spectra(i).Title=path{i};

% Clear temporary variables
clearvars -except path spectra
end
%plot
scatter(spectra(1).ppm1,spectra(1).ppm2,4)
