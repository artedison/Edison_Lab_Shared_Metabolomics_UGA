%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     WORKFLOW for 1D HRMAS Samples
%     MTJ 2021

%{

    The raw data can be found here: 
            > 
    The data were processed in NMRPipe, following the
        using the files in:
            > nmrPipe processing for HRMAS_ananaerobic_13c_1 - Notes.txt
        which is found in the processed data folder. 
    
    All data and results and metadata from this workflow will be stored here:
            > 
        and in 
            > 
        and in
            > 
    
    Original NOESY data on NMR Synology: /NMRbackup/carbon600/nmrdata/judgemt/HRMAS_ananaerobic_13c_1
    
    
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd(scripts_directory);save([sampleName,'_noesy.mat'])
% This command will load the current saved .mat file for the workflow:

%load([sampleName,'.mat')

%% STEP 1: NAVIGATE TO THE FOLDER IN WHICH THIS FILE SITS
    % Also, add the Edison Lab Metabolomics Toolbox:
%     addpath(genpath('/Users/mjudge/Edison_lab_UGA'))

%%  
    cd('mFileLocation')
    expName.sampleName = 'sampleNameGoesHere';
    
    [~,~,thisFile] = findCurrentFile();
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    %% Load the Processed data from ft files
        % Load all the ft files as spectra
            cd ftDataPath_relative
                loadallft();
            expName.spectra = spectra;
            
        % Get timepoints from acqu files
            expName.dataDir = {spectra.Title}';
            cd rawDataPath_relative
                    [expName.startTimes,expName.dateTimes] = getRunTimes_NMR(expName.dataDir);
                    expName.timePoints = expName.startTimes - expName.startTimes(1); % these are in hours
            cd(project_directory) 
         
%         figure, plot(expName.ppm,expName.X)
%             set(gca,'XDir','reverse')
%             xlabel('Chemical Shift (ppm)')
%             ylabel('Signal Intensity')

    %% Referencing (if you have a calibration peak) 
        
        expName.spectraRef = expName.spectra; % in case you don't want to reference
        % pyruvate @ 2.364ppm
        % peak @ 1.6ppm
        expName.ref_threshold = .009;
        expName.ref_offset = 1.6;
        expName.spectraRef = ref_spectra(expName.spectra,expName.ref_threshold,expName.ref_offset);
        %expName.spectraRef = ref_spectra(expName.spectra,expName.ref_threshold,expName.ref_offset,'testThreshold');
        close(gcf)
        expName.spectraRef = expName.spectra;
        
    %% Make 2D matrix of 1D spectra using sorted spectra

        [expName.X,expName.ppm,~]=Setup1D(expName.spectraRef); 
            %expName.ppm = expName.ppm + 2.364;
            
    %% Chop off at n hours
    %     cutTime = 13;
    %     [~,cutInd] = min(abs(times_13c1d - cutTime));
    %     X_13c1d = X_13c1d(1:cutInd,:);
    %     times_13c1d = times_13c1d(1:cutInd);

    %% Remove ends, water region
    
        expName.X = remove_region(expName.X,expName.ppm,4.7,5);
        [expName.X,expName.ppm] = remove_ends(expName.X,expName.ppm,-0.5,10);
        
%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %% Normalize (to DSS, if applicable)  
    
%         matrix = expName.X;
%         currentppm = expName.ppm;
%         expName.NormMethod = 'maximum'; % alternative, 'integral'
%         expName.ROInorm = []; % empty indicates that the peak region will be selected 
% 
%         [expName.XN,expName.ROInorm,~] = normalize_HRMASdata(matrix,currentppm,expName.NormMethod,expName.ROInorm);

     %% Collapse spectra        
        if isfield(expName,'XN')
            matrix = expName.XN;
        else
            matrix = expName.X;
        end   
        expName.binsize = 10;
        %[Xcollapsed,timesCollapsed,totalTime,resolution] = HRMAS_collapseTimes(matrix,timePoints,binsize);
        [expName.Xcollapsed,expName.timesCollapsed,expName.totalTime,expName.resolution] = HRMAS_signalAverage(matrix,expName.timePoints,expName.binsize);
    
     %% Make a Stack Plot of the spectra:
        matrix = expName.Xcollapsed(1:10:end,:);
        %matrix = expName.X;
        currentppm = expName.ppm;
        
        expName.horzshift = -0.007;%0;
        expName.vertshift = 10E7;
        expName.xlims = [0.868577188940092,1.31500576036866];
        expName.ylims = [-1866365131.57895,452713815.789474];

        expName.plotTitle = {['Time-averaged ','expName',' data for ',expName.sampleName],[expName.resolution,'-min resolution, ',expName.totalTime,'-h total']};

        stackSpectra(matrix,currentppm,expName.horzshift,expName.vertshift,expName.plotTitle)
        set(gcf, 'InvertHardCopy', 'off');
%             set(gca,'xlim',expName.xlims)
%             set(gca,'ylim',expName.ylims)
        %xlim([0.5,5.2]) % in ppm, for auto-zoom to a region of interest
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    