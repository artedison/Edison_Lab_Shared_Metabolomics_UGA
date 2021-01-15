function [output] = constructHRMASDirectory_4(destinationDir,newDataDir)
%%
% Note: replaced by constructHRMASDirectory.m

% Next: don't overwrite existing files


%% ***** Rework to conform to the sample(s).expType(t) format

%% ***** Rework paths storage (all under paths)

%%
% Make sure destination exists

    if ~exist(destinationDir,'dir')
        mkdir(destinationDir);
    end

% Get Sample Names and Locations

      newDataDirs=dir(newDataDir);
      zipFiles = newDataDirs(contains({newDataDirs.name},{'.tar.gz','.zip'}));
      newDataDirs = newDataDirs(~contains({newDataDirs.name},{'.','..','.DS_Store','.tar.gz'}));
      
% Check to make sure the newDataDir isn't empty
      if isempty(newDataDirs)
          warning(['The data input directory, ',newDataDir,', is empty. Please copy your datasets to there and re-run. Aborting constructHRMASDirectory_4().'])
          return
      end
      
% For each sample:

    sample = struct();
    
    for s = 1:length(newDataDirs)
        
    % Get the Sample Name
    
        sample(s).name = newDataDirs(s).name;     
        
    % Make the standard new directories and name them
            
            cd(destinationDir)
                mkdir(sample(s).name)
                
                cd(sample(s).name)
                    sample(s).paths.sample = pwd;
            
                    mkdir('data');mkdir('results');mkdir('scripts');
                        sample(s).paths.data = [pwd,'/data'];
                        sample(s).paths.results = [pwd,'/results'];
                        sample(s).paths.scripts = [pwd,'/scripts'];
                        
                    cd('data');mkdir('raw');mkdir('processed');  
                        sample(s).paths.raw = [pwd,'/raw'];
                        sample(s).paths.processed = [pwd,'/processed'];
                                cd('./processed');mkdir('nmrpipe');
                                sample(s).paths.nmrpipe = [pwd,'/nmrpipe'];
            sample(s).paths.originalData = [newDataDir,'/',sample(s).name];

%             cd(sample(s).paths.scripts)
%                     mkdir('nmrPipe_templates')
%                     cd('nmrPipe_templates')
%                         sample(s).paths.templates = pwd;
%                         mkdir('representative_spectrum')                        

    % Sort the Experiments into Separate Directories by Experiment Type
        
        
        % Generate map: <expType types> - filenames
        
            % Get the acqus files (in the order they are read; NOT natural number sorted)
                
                sample(s).paramFiles.name = sample(s).name;
                sample(s).paramFiles.fileData = readAcqusFiles(sample(s).paths.originalData);
                
        % Define the <expType types> 
               % Find types of expTypes from that list
                    sample(s).expTypes = unique({sample(s).paramFiles.fileData.experimentType});
                    
                    sample(s).dataFiles = dir(sample(s).paths.originalData);
                        sample(s).dataFiles(ismember({sample(s).dataFiles.name},{'.','..','.DS_Store'})) = []; % remove unnecessary dirs                             
                        
        % Before doing anything else, check to make sure the data exist:
            
            % Index for empty files
                neededFiles = {'acqus','fid','pulseprogram','pdata'};
                neededpdata = {'title'}; % 1i, 1r could be added. sometimes these don't make it in
                neededAll = [neededFiles,neededpdata];
                sample(s).requiredFiles = neededAll;
                checklist = zeros(length(sample(s).dataFiles),length(neededAll));
                
                for d = 1:length(sample(s).dataFiles)
                    sample(s).dataFiles(d).fullname = [sample(s).dataFiles(d).folder,'/',sample(s).dataFiles(d).name];
                    expFiles = dir(sample(s).dataFiles(d).fullname);
                    pdataFiles = dir([sample(s).dataFiles(d).fullname,'/pdata/*/']);
%                                         expFiles = dir(fullfile(sample(s).dataFiles(d).fullname,'**/*.*'));
                    checklist(d,:) = [ismember(neededFiles,{expFiles.name}),...
                                             ismember(neededpdata,{pdataFiles.name}),...
                                             ];
                    sample(s).dataFiles(d).hasData = all(checklist(d,:));
                    sample(s).dataFiles(d).missingData = neededAll(   ~checklist(d,:)   );

                end
            % Store info about missing data
                sample(s).filesMissingData.paramFiles = sample(s).paramFiles.fileData(~[sample(s).dataFiles.hasData]);
                sample(s).filesMissingData.dataFiles = sample(s).dataFiles(~[sample(s).dataFiles.hasData]);
                sample(s).filesMissingData.missing = {sample(s).dataFiles(~[sample(s).dataFiles.hasData]).missingData}';
                
            % Clear the paramFiles and dataFiles lines for empty data files
                sample(s).paramFiles.fileData(~[sample(s).dataFiles.hasData]) = [];
                
                sample(s).dataFiles(~[sample(s).dataFiles.hasData]) = [];

            
            % Alert the user that datapoint(s) will be missing, if
            % necessary
                if ~isempty(sample(s).filesMissingData.missing)
                    warning(['Some files in   ',sample(s).name,'   did not contain necessary data. See info in studyInfo.sample(',num2str(s),').filesMissingData.'])
                end

            % Clear the fields without data
            
                    sample(s).dataFiles = rmfield(sample(s).dataFiles,'missingData');
                    sample(s).dataFiles = rmfield(sample(s).dataFiles,'hasData');

        % Each expt needs its own templates
            cd(sample(s).paths.scripts)
                for t = 1:length(sample(s).expTypes)
                    sample(s).expType(t).type = sample(s).expTypes{t};
                    sample(s).expType(t).paths.scripts = pwd;
                        mkdir(sample(s).expType(t).type), cd(sample(s).expType(t).type)   % MTJ edit 21DEC2020 to allow different processing dir for each expt.
                        mkdir('nmrPipe_templates')
                        cd('nmrPipe_templates')
                            sample(s).expType(t).paths.templates = pwd;
                            mkdir('representative_spectrum') % no need to add to paths list
                end
        % Make a directory for each <expType type> within "raw"
            cd(sample(s).paths.raw)
                for t = 1:length(sample(s).expTypes)
                    % Make the directory for this type
                        
                        mkdir(sample(s).expType(t).type)
                        sample(s).expType(t).paths.raw = [pwd,'/',sample(s).expType(t).type];
                end       
                
        % Make fid and ft directories for each <expType type> within "nmrpipe"
                for t = 1:length(sample(s).expTypes)
                    
                    cd(sample(s).paths.nmrpipe)
                        mkdir(sample(s).expType(t).type)
                        cd(sample(s).expType(t).type)
                            mkdir('fid'); sample(s).expType(t).paths.fid = [pwd,'/fid'];
                            mkdir('ft');  sample(s).expType(t).paths.ft = [pwd,'/ft'];
                            mkdir('proc_files'); sample(s).expType(t).paths.proc_files = [pwd,'/proc_files'];
                            cd('proc_files');
                                mkdir('fid_com'); sample(s).expType(t).paths.fid_com = [pwd,'/fid_com'];
                                mkdir('ft_com');  sample(s).expType(t).paths.ft_com = [pwd,'/ft_com'];      
                                
                end  
                cd(sample(s).paths.nmrpipe)
                %sample(s).expType(t) = sample(s).expType; 
                
        % Copy (Move) the expTypes and zip/tar file for each <expType type>
        
                % Index the spectra based on their expType
                
                    [~,sample(s).exptKey] = ismember({sample(s).paramFiles.fileData.experimentType},sample(s).expTypes);
        
                for t = 1:length(sample(s).expTypes)
                    
                    % Find the file indices for this type
                    
                        inds = sample(s).exptKey == t;
                        
                    % Use system commands to move the files of that type
                    
                        cd(sample(s).paths.originalData); % move from original data location
                        
                            % use escape characters for dropbox path
                                cmdLineExpdir = regexprep(sample(s).paths.sample,' ','\\ ');
                                cmdLineExpdir = regexprep(cmdLineExpdir,'(','\\(');
                                cmdLineExpdir = [regexprep(cmdLineExpdir,')','\\)'),'/data/raw/',sample(s).expType(t).type];
                            %cmd = ['cp -R ', sprintf('%s',sample(s).dataFiles(inds).name), cmdLineExpdir]; % doesn't work/too slow
                                                            
                        cmd = ['mv ', sprintf('%s ',sample(s).dataFiles(inds).name), cmdLineExpdir];
                        system(cmd);
                        
                        sample(s).expType(t).files = sample(s).dataFiles(inds);
                        
                end  
                
                % Clean up old dir
                    cd([sample(s).paths.originalData]),cd ..
                    system(['rm -r ',sample(s).name])
                    
                cd(sample(s).paths.raw) % Go back to the raw data dir 
                
            fprintf(['\n\tSample ''',sample(s).name,''' constructed successfully...\n\n'])
    end
    
    output.zipFiles = zipFiles;
    output.sample = sample;
    cd(sample(1).paths.sample),cd ..
    fprintf('\n\tconstructHRMASDirectory_4() completed. Revelant data are found in the output struct. Run processCIVMdata() next.\n\n')
end

