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
for i=1:length(path)
[~, ~, raw] = xlsread('/Users/chaevienclendinen/Dropbox/Chaevien/13C__elegans/INAD/lightintense.xlsx',path{i});
spec = raw(2:end,1:3);
% spec(cellfun(@isnan, spec))=[];
% spec=spec';
data(:,:) = reshape([spec{:}],size(spec));

% array to column variable names
ppm1 = data(:,2);
DQ = data(:,3);
Inten = data (:,1);

% Make Spectra structure

spectra(i).ppm1=ppm1;
spectra(i).ppm2=DQ;
Intensity=zeros(400,400);
Intensity(sub2ind(size(Intensity),round(ppm1), round(DQ)))=data(:,1);
spectra(i).real=Intensity;
spectra(i).Title=path{i};

% Clear temporary variables
clearvars -except path Intense
end

%% plot
scatter(spectra(1).ppm1,spectra(1).ppm2,6)
set(gca,'XDir','reverse','YDir','reverse')

%%

    
    
    