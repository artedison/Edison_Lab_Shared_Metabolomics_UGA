function [data] = HRMAS_nmr_runStdProc(sampleInfo,sample,expType)
%% HRMAS_nmr_runStdProc
% 
% Run standard basic processing on a dataset (easy to put in a loop). Facilitates
% concatenation (hopefully). Built from HRMAS_nmr_processAndStoreData4.m
% 
% MJ 13APR2018
% YW 10/04/2018 add local storing 
% MJ 10/2020 rework for processCIVMdata()

%% handle inputs

specList = sampleInfo.sample(sample).expType(expType);

%% Take care of filepath stuff

    data_directory = specList.paths.ft;
    rawData = specList.paths.raw;
    
%% PROCESSING
    %% Load the Processed data from ft files
        % Load all the ft files as spectra
            cd(data_directory)
            loadallft();
        % Get timepoints from acqu files
            dataDir = {spectra.Title}';
            cd(rawData)
            [startTimes,~] = getRunTimes_NMR(dataDir);
            timePoints = startTimes - startTimes(1); % these are in hours
        % Get pulsewidths too
            [P0] = getP0_NMR(dataDir);            

    %% Referencing
        spectraRef = ref_spectra(spectra,.004,0);
        close(gcf)

%% STORAGE

    data.spectra = spectraRef;  
    % Add in fields for sample and experiment so we can concatenate with
    % other data later on
        [data.spectra(:).sample] = deal(sampleInfo.sample(sample).name);
        [data.spectra(:).experiment] = deal(sampleInfo.sample(sample).expTypes{expType});
    for i = 1:length(spectra)
        data.spectra(i).startTime = startTimes(i);
    end
    data.pulseWidthP0 = P0;
    data.timePoints = timePoints;
    
end
