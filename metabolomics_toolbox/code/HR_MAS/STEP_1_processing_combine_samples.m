%% Workflow for combining Matlab-processed samples in an HR-MAS run
% This workflow assumes no normalization or collapsing.
%
% NMRPipe Processing:
% - conversion to .fid
% - phasing
% - poly baseline correction
% - fourier transform (conversion to .ft)
% Matlab Processing:
% - batch referencing to DSS
% - remove water, ends
%     - really want to do this? ppms might be different?
%     - this is done in this script. Not included in the shorter ones
% 'HRMAS_ncrassa_paper_Sample_[samples]_1h1d_short.m'

% Before running, add the following to the Matlab path:
% path of datafolder (analysis folder)

% The following is an example of how to add the Edison Lab shared toolbox (Edison_Lab_Shared_Metabolomics_UGA) to the
% local path (REQUIRED).
    % addpath(genpath('/Users/mjudge/Edison_lab_UGA'))

% All samples were run in the same media, same wash media, same growth
% media. Same conidia over the course of ~3wks

%% Orient in Directory
toolboxfolder='/Users/yuewu/Documents/GitHub/Edison_Lab_Shared_Metabolomics_UGA/';
addpath(genpath(toolboxfolder));%%the path for toolbox
anafolder='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/archive/test.code/analysis/';
cd(anafolder);
addpath(genpath(anafolder)); %% the path for datafolder
thisFile='STEP_1_processing_combine_samples.m';
wd=which(thisFile);cd(wd(1:end-length(thisFile)-1));clear('wd');
meta_directory=cd();
% mkdir('results')
cd 'results', results_directory=cd();

%% Reprocess data manually (if necessary)
samples=[4,5,6,10,11,12]; %here only 6 replicates are processed for example, while in the dataset we provided, there are 9 replicates.
sampleKey={'aerobic',	'aerobic',	'aerobic',	'oxygen-limited',	'oxygen-limited',	'oxygen-limited'};
%Open all processing files in Matlab (for editing)
% for s = samples
%     eval(sprintf('open(''HRMAS_ncrassa_paper_Sample_%d_1h1d_short.m'')',s))
% end

% run the scripts under each data folder
% for s = samples(1:6)
%         eval(sprintf('HRMAS_ncrassa_paper_Sample_%d_1h1d_short',s))
%         fprintf(['\n\tWriting the .mat file for Sample ',num2str(s),' ...\n\n']);
%         cd(meta_directory);save(['HRMAS_ncrassa_paper_Sample_',num2str(s),'_data.mat']);
%         clear('data_directory','dataDir','dateTimes','P0','ppm_1h1d','project_directory','scripts_directory','spectra','spectraRef_1h1d','startTimes','thisFile','times_1h1d','X_1h1d')
% end

%% Process all the data using HRMAS_nmr_processAndStoreData
% This is the recommended method.
% modify HRMAS_nmr_processAndStoreData.m as necessary (applies params to
% all samples)

% future edits will include nmrPipe processing
filenameTemplate='HRMAS_ncrassa_paper_Sample_';
for s=samples
    HRMAS_nmr_processAndStoreData2([filenameTemplate,num2str(s)],'noSave',true);
    fprintf(['\n\tWriting the .mat file for Sample ',num2str(s),' ...\n\n']);
    cd(meta_directory);
    % save([filenameTemplate,num2str(s),'_data.mat']);
    clear('thisFile')
end
%% Compare collapsing in Matlab with fid-averaging in nmrPipe
%  HRMAS_ncrassa_paper_testingTimeAveraging.m
    % (The results are the same.)

%% Get the Matlab-processed data:
sampleData=struct();
for s=1:length(samples)
    eval(sprintf('load(''HRMAS_ncrassa_paper_Sample_%d_1h1d_data.mat'',''X_1h1d'',''ppm_1h1d'',''times_1h1d'')',samples(s)))  % aerobic
    sampleData(s).X_1h1d=X_1h1d;
    sampleData(s).ppm_1h1d=ppm_1h1d;
    sampleData(s).times_1h1d=sort(times_1h1d);
    sampleData(s).sampleNumber=samples(s);
    sampleData(s).condition=sampleKey{s};
end
clear 'ppm_1h1d' 's' 'times_1h1d' 'X_1h1d' 'sampleKey'

%% Further Matlab Processing in batch

%% Remove ends, water region
    % Generate an average
        % Ends
ends=[-0.5,10];
for s=1:length(samples)
    [sampleData(s).XR_1h1d, sampleData(s).ppmR_1h1d]=remove_ends(    sampleData(s).X_1h1d,sampleData(s).ppm_1h1d, ends(1),ends(2));
end
        % Water
waterBounds=[4.7,5];
for s=1:length(samples)
    sampleData(s).XR_1h1d=remove_region(    sampleData(s).XR_1h1d,sampleData(s).ppmR_1h1d,  waterBounds(1),waterBounds(2));
end
%% Chop off at n hours (optional)
cutTime=11;
for s=1:length(samples)
    [~,cutInd]=min(abs(sampleData(s).times_1h1d- cutTime));
    sampleData(s).XR_1h1d= sampleData(s).XR_1h1d(1:cutInd,:);
    sampleData(s).timesR_1h1d= sampleData(s).times_1h1d(1:cutInd);
end

%% Normalize (to DSS)
ROInorm=[-0.5,0.5]; % DSS bounds
method='maximum'; % use intensity at DSS peak maximum as norm. factor
for s=1:length(samples)
    [sampleData(s).XN_1h1d,sampleData(s).DSS_ROI,  sampleData(s).normFactors]=normalize_HRMASdata(   sampleData(s).XR_1h1d,sampleData(s).ppmR_1h1d,    method,ROInorm);
end

%% Collapse spectra
binSize=3;
for s=1:length(samples)
    [sampleData(s).Xcollapsed_1h1d,sampleData(s).timesCollapsed_1h1d,sampleData(s).totalTimeCollapsed,sampleData(s).resolutionCollapsed]=       HRMAS_collapseTimes(sampleData(s).XN_1h1d, sampleData(s).timesR_1h1d,binSize);
end

%% Go back and get the pulse widths for each run:
% Add all the directories to the path
for i=1:length(samples)
    % Find the file
    eval(sprintf('thisFile = ''HRMAS_ncrassa_paper_Sample_%d_1h1d_short.m'';',samples(i)));
    wd=which(thisFile);cd(wd(1:end-length(thisFile)-1));clear('wd')
    cd ..
    eval(sprintf('dataDir = ''data/NMR/HRMAS/HRMAS_ncrassa_paper_Sample_%d/HRMAS_ncrassa_paper_Sample_%d'';',samples(i),samples(i)))
    % Get p0 from acqu files
    sampleData(i).P0=getP0_NMR(dataDir);
end
cd(meta_directory)
clearvars -except red yellow meta_directory results_directory sampleData samples

%% Save the Workspace
save('sampleData.mat','sampleData','samples','meta_directory','results_directory')

%% Tracing Ridges
% See: STEP_2_ridge_tracing.m
