%% iSTOCSY_mj
% % 2FEB2017
% Get to the starting directory
% Make sure you have an 'Inputs' and an 'Outputs' folder!
% 
% startDir:     string format address of working directory in which outputs
%               are desired
% steps:        number of rounds (iterations) of STOCSY, where the outputs
%               (drivers and responders) for the first STOCSY are fed back into 
%               the 
% X:            spectral matrix (full spectral matrix is preferred for maximum
%               variance
% ppm:          ppm matrix, from low to high ppm values
% targets:      ppm values to be STOCSY'd (also called 'drivers')
% threshold:    STOCSY correlation threshold magnitude (0-1) for determining interactions.
%               STOCSY_mj_copy handles positive/negative correlations
% intT:         threshold for determining minimum peak height for counting correlation points as interactions 
% figureOption: 'generateFigures' pops out, saves (.fig), and closes a
%               customized figure, showing 1) Raw STOCSY output 2) Thresholded STOCSY 
%               3) All spectra with axes linked for easy investigation, for each driver in the run. 
%               This takes a lot more time and memory. Entering anything else, like 'noFigures', or
%               'enerateFigures' will suppress this and run MUCH faster.
%               Suppression recommended for runs when determining
%               thresholds, or just getting interactions. The function can
%               always be run on smaller lists of more interesting peaks as
%               primary drivers, if STOCSY figures are desired.
% setup:        'peakpick' or 'peaklist'. If 'peakpick', then the
%               Peakpick1D function is run on the data to generate a list 
%               of drivers. This is desirable if mapping and visualizing 
%               correlation structure (recommended), and parameters must be
%               adjusted within the call to Peakpick1D(). If there is a list of 
%               interesting peaks (or one), use 'peaklist' but be aware
%               that an input file (.csv) is required, and must be provided
%               for      
%               peaklist: conditional variable (not necessary under
%                           'peakpick', but necessary under
%                           'peaklist'). This should be a vector containing
%                           the ppms of the peaks for which STOCSY is
%                           desired.
%{
 Outputs and inputs are stored in directories under "iSTOCSY" in startDir,
 the file structure typically looks like this:
 
iSTOCSY: (this is the main directory
    ->Inputs (self-explanatory)
        'primaryDrivers_picked_peaks_' [timestamp] '.csv'
            ^this file will be made programmatically^
        or
        'primaryDrivers_customInput_' [timestamp] '.csv'
            ^this file will be made programmatically^
            and
        '[peakList].csv'
            ^this file will need to be supplied manually if 'peakpick' option
            is not used.
    ->Outputs:

        ->Steps_2_Threshold_0.9_[timestamp] 
%{
                this is a subdirectory made for this particular instance (run) 
                of the program. "Steps_2" indicates a two-step iSTOCSY,
                "Threshold_0.9" indicates a STOCSY cutoff of +/- 0.9
                correlation coefficient. [timestamp] is a date-time number
                to make sure directories are trackable and unique (no
                overwriting old data).
 %}
                                    
            ->STEP_1_threshold_0.9

                STOCSY_clusters_STEP_1_threshold_0.9.csv
%{
                        this file can be opened and individual rows contain
                        the ppms that are the max of each ppm range that
                        was correlated with the driver, which is in the
                        first column. One-to-many relationship. Can be used
                        in COLMAR 1H trace searches against BMRB,HMDB,etc.
                        databases. Best to overlay results with your own
                        data for verification by visual matching.
%}

                STOCSY_clusters_STEP_1_threshold_0.9_binaryInteractions.csv
%{
                        this file can be directly uploaded into Cytoscape for
                        visualization, one to one relationships.
                        Correlations in third column.
                            - New session -> With Empty Network -> [name]
                            -> File -> Import -> Network -> File [select file, open]
                            - Click Column 1 header -> Data Type->string ->
                            Meaning->Source Node
                            - Click Column 2 header -> Data Type->string ->
                            Meaning->Target Node -> okay => Network
                            - Click Column 3 header -> Data Type->decimal
                            -> Meaning-> Edge Attribute -> okay => Network
                        To be useful in representing correlation structure, edges 
                        may additionally be colored/shortened using correlation 
                        coefficient. Highly connected clusters will be more reliable
                        for dereplication. Star topologies are usually not
                        reliable, but you must run >1 steps to avoid these.
                            - MCODE with haircut
                            - I like to change the styles and save a
                            preferred one. Dark background easier to look
                            at. 
                            - To get rid of edge labels: Style tab -> 
                            - My preferred layouts:
                                - Layout -> yFiles -> organic (shows
                                connectivity/topology type well
                                - Layout -> Prefuse Force-Directed Layout
                                (pretty, looks biological)
                                - Layout -> Edge-weighted Spring Embedded
                                Layout -> select the correlations (computationally intensive, will
                                take 30min-1h for large networks (~90K
                                interactions). Shows Connectivity and
                                stronger correlations as more compact.
%}
                Figures for STEP_1 will be saved to here as well, if
                activated by 'generateFigures'.

                ->STEP_2_threshold_0.9

                    STOCSY_clusters_STEP_2_threshold_0.9.csv

                    STOCSY_clusters_STEP_2_threshold_0.9_binaryInteractions.csv
                    
                    Figures for STEP_2 will be saved to here as well, if
                    activated by 'generateFigures'.

                    -> STEP_3....
%}


% {
%% Setup (for testing):
    load WorflowRMB_1D_1DSTOCSY_MJ_15JUN2016
    [~,~,raw] = xlsread('Decoder_key_JG-to-RBMJ_reorder_for_replicates.xlsx');
    X = XRALN([1:2,4:8,10:38],:); %Get rid of media blanks, pooled samples (> 39), weird samples (28:30)

    clearvars('-except','X','XRALN','ppm');
%}
%% Input variables
    startDir = '/Users/mjudge/Documents/MATLAB/Edison_lab_UGA/MatLab_RMBorges/Not_On_RMBMac';
    steps = 1;
    X = XRALN;
    threshold = 0.90;
    intT = 2E-4; % intensity threshold; if responding 'peaks' are less intense than this in the mean of all X spectra, then set to zero (cut out picking up noise; get this manually)
    figureOption = 'generateFigures'; %
        %NOTES: generating figures will also save them in the outputs folder. Type anything else, like 'noFigures' to run fast and just get the interactions.
        %       targets = []; % or, just enter from within Matlab. Better to use
        %       files as records, though.
    %ppm = ppm;
%% Choose setup
setup = 'peakpick'; 
%setup = peaklist; peaklist = ''.csv;

switch setup
    case 'peakpick'
        
        %% STANDARD SETUP (PEAKPICKING)
            % Set your starting (working) directory. This should contain an 'Inputs'
            % folder, where you keep the peaklist you want to use.
                cd(startDir)
                mkdir(startDir,'iSTOCSY')
                cd([startDir '/iSTOCSY'])
                mkdir('iSTOCSY_Inputs') 
                mkdir('iSTOCSY_Outputs')
                % use this to make the inputs directory
                % A warning will appear after this command if the directory already exists, but a pre-existing directory
                % with that name will not be overwritten. One way of checking for existence.
                inDir = [startDir '/iSTOCSY' '/iSTOCSY_Inputs'];
                outDir = [startDir '/iSTOCSY' '/iSTOCSY_Outputs'];
                addpath(genpath([startDir]))
                cd(inDir)
                
                % Do Peakpicking for drivers. We want good coverage of the spectrum.
                    [pp.peakmatrix,pp.shifts]=Peakpick1D(X,ppm,'mean',.7,'Complex');
                        %[pp.peakmatrix,pp.shifts]=Peakpick1D(X,ppm,'max',.5,'Complex'); 
                        %maybe save the figure for future reference?
                
                % Write the results to file for safe keeping.
                primaryDrivers = ['peak_picked_X_' num2str(now) '.csv'];
                fprintf(['Creating \"' primaryDrivers '\" as primary drivers file in "Inputs" directory...\n'])
                csvwrite(primaryDrivers,pp.shifts');
                fprintf(['Using \"' primaryDrivers '\" as primary drivers file...\n'])

            % Navigate to Starting directory 
                cd([startDir '/iSTOCSY'])

    case 'peaklist'
        %% SETUP FOR PRE-EXISTING DRIVER PEAK LIST - NOT TESTED
            % Do Peakpicking for drivers. We want good coverage of the spectrum.
                cd(startDir)
                mkdir(startDir,'iSTOCSY')
                cd([startDir '/iSTOCSY'])
                mkdir('iSTOCSY_Inputs') 
                mkdir('iSTOCSY_Outputs')
                % use this to make the inputs directory
                % A warning will appear after this command if the directory already exists, but a pre-existing directory
                % with that name will not be overwritten. One way of checking for existence.
                inDir = [startDir '/iSTOCSY' '/iSTOCSY_Inputs'];
                outDir = [startDir '/iSTOCSY' '/iSTOCSY_Outputs'];
                addpath(genpath([startDir '/iSTOCSY']))
                cd(inDir)

                primaryDrivers = csvread(peaklist);

                fprintf(['Creating \"' primaryDrivers '\" as primary drivers file in "Inputs" directory...\n'])
                csvwrite(primaryDrivers,pp.shifts');
                fprintf(['Using \"' primaryDrivers '\" as primary drivers file...\n'])

            % Navigate to Starting directory 
                cd(startDir)
end
%% You need the following in order to run in most situations:
    % Check all variables before running
%% iSTOCSY_mj variables      

% Start with a list of drivers, Then run STOCSY:
    drivers = csvread(primaryDrivers);

% Make output directory
    cd(outDir)
    STOCSYoutputs = [ 'Steps_' num2str(steps) '_Threshold_' num2str(threshold) '_' num2str(now)];
    mkdir(STOCSYoutputs);
    cd(STOCSYoutputs);

    interactions = drivers(1:3); % just for the first run

%% Run the iSTOCSY_mj loop
    for currentStep = 1:steps
        % Set Up for this run
            % Convert interactions to list of drivers
           drivers = interactions2drivers(interactions);

            % make the directory for Outputs, cd there. Effect is nested folder
            % structure.
                outputFolder = ['STEP_' num2str(currentStep) '_threshold_' num2str(threshold) ];
                mkdir(outputFolder);
                cd(outputFolder);
                str = ['STEP_' num2str(currentStep)];

        % Run all the STOCSYs
            eval(   sprintf('STOCSY_output_step_%d = STOCSY_mj_wrapper(X,ppm,drivers,threshold,intT,figureOption,str);', currentStep)   );
                % this little line dynamically creates the STOCSY_output_step_
                % structure so that the step number is recorded.
            % If 'interactions field in the results structure' is used, read in from file, zeros get in
            % because Matlab can't handle uneven rows and fills with zeros.
            % Instead, take the same numbers from BinaryInteractions(:,1:2), which
            % has no zeros. This was a real bugger because it threw errors
            % in deeper functions, while the issue was up here! 
 
            eval(   sprintf('interactions = STOCSY_output_step_%d.BinaryInteractions(:,1:2);', currentStep)   );
            interactions = unique(sort(reshape(interactions,[],1)));
        %}
            
    end
clear('');
cd(startDir)
save('workspace')
%% 

    
    
    
    
    
    
    
    
    
    
    