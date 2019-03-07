function HRMAS_nmr_processAndStoreData4(sampleName,fileName,expType,saveData,local)
%dirName = 'HRMAS_ncrassa_paper_Sample_4'
%saveData = 'save';
%local=false;

% MJ 13APR2018
% YW 10/04/2018 add local storing 

%%
    % Assume in parent directory with all samples 
    %localdir=cd(fileName);

%% Take care of filepath stuff
    cd(sampleName)       
    project_directory = cd();
    cd scripts
        scripts_directory = cd();
        cd(project_directory)
    cd data
        data_directory = cd();
        cd(project_directory)

%% PROCESSING
    %% Load the Processed data from ft files
        % Load all the ft files as spectra
            cd(data_directory)
            cd(['nmrpipe/',expType,'/ft'])
            loadallft();
        % Get timepoints from acqu files
            dataDir = {spectra.Title}';
            cd(['../../../raw/',expType])
            [startTimes,~] = getRunTimes_NMR(dataDir);
            [P0] = getP0_NMR(dataDir);
            timePoints = startTimes - startTimes(1); % these are in hours
            cd(project_directory); cd results

    %% Referencing
        spectraRef = ref_spectra(spectra,.004);
        close(gcf)

    %% Make 2D matrix of 1D spectra using sorted spectra
        % Get the sorted spectra in a format that Setup1D understands
            [X,ppm,~]=Setup1D(spectraRef); % see that X rows are in the order of the natural-number-sorted filenames.
            % Remove the last timepoint, as this is usually incomplete:
                %X = X(1:end-1,:);
                %timePoints = timePoints(1:end-1);

    %% Save
%         if local
%           scripts_directory=localdir;
%         end
        if strcmp(saveData,'save') || local
            fprintf(['\n\tWriting "',fileName,'_data.mat" in \n\t',scripts_directory,'...\n\n'])
            cd(scripts_directory);save([expType,'_data.mat'])
        end
end
