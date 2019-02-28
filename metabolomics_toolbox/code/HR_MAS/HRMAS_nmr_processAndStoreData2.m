function HRMAS_nmr_processAndStoreData3(dirName,saveData,local)
%dirName = 'HRMAS_ncrassa_paper_Sample_4'
%saveData = 'save';
%local=false;

% MJ 13APR2018
% YW 10/04/2018 add local storing 

%%
    localdir=cd();
    thisFile = [dirName,'_1h1d_short.m']
    wd = which(thisFile);cd(wd(1:end-length(thisFile)-1));clear('wd')

%% Take care of filepath stuff
    cd ..        % go up one directory
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
            cd(['NMR/HRMAS/',dirName,'/nmrpipe/ft_1h1d'])
            loadallft();
        % Get timepoints from acqu files
            dataDir = {spectra.Title}';
            cd(['../../',dirName])
            [startTimes,dateTimes] = getRunTimes_NMR(dataDir);
            [P0] = getP0_NMR(dataDir);
            times_1h1d = startTimes - startTimes(1); % these are in hours
            cd(project_directory); cd results/1h1d

    %% Referencing
        spectraRef_1h1d = ref_spectra(spectra,.004);
        close(gcf)

    %% Make 2D matrix of 1D spectra using sorted spectra
        % Get the sorted spectra in a format that Setup1D understands
            [X_1h1d,ppm_1h1d,~]=Setup1D(spectraRef_1h1d); % see that X rows are in the order of the natural-number-sorted filenames.
            % Remove the last timepoint, as this is usually incomplete:
                X_1h1d = X_1h1d(1:end-1,:);
                times_1h1d = times_1h1d(1:end-1);

    %% Save
        if local
          scripts_directory=localdir;
        end
        if strcmp(saveData,'save') || local
            fprintf(['\n\tWriting "',dirName,'_1h1d_data.mat" in \n\t',scripts_directory,'...\n\n'])
            cd(scripts_directory);save([dirName,'_1h1d_data.mat'])
        end
end
