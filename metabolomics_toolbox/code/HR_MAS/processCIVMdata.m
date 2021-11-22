function [updatedStudyInfo] = processCIVMdata(studyInfo,destinationDir,newDataDir,templateDir,varargin)
%% processCIVMdata(studyInfo,varargin)
%
%     This function provides an interactive way to process single or CIVM datasets
%     A CIVM dataset is defined as a time series run collected on a distinct sample,
%     Examples of a sample include the following standard CIVM controls:
%         - Culture Media (~10 sequential spectra)
%         - Culture Media + Spike-in carbon source  (~10 sequential spectra)
%         - Culture Media + Organism (~1 h of sequential spectra, 100's of spectra)
%     as well as the "actual" run:
%         - Culture Media + Organism + Spiked carbon source (hours of collection, 1000's of spectra)
%
%     When a dataset is recorded in TopSpin (via iconNMR or multizg, etc.), the
%     spectra are recorded under a uniquely named directory. For instance:
%         "civm_ncrassa_qax_03_media"
%
%     For the experiment structure listed above, the directories would be:
%         civm_ncrassa_qax_03_media
%         civm_ncrassa_qax_03_media_org
%         civm_ncrassa_qax_03_media_qa
%         civm_ncrassa_qax_03_media_org_qa_iconRun
%
%     These will all be kept under civm_ncrassa_qax_03, which is the name
%     of this overall experiment.
%
%     Individual spectrum directories are contained under each "sample directory"
%     or dataset. One of these files, i.e. "Bruker Files", for timepoint 1
%     would look like this:
%
%             1
%                 acqu
%                 acqus
%                 ...
%                 pdata
%                     1
%                         1i
%                         1r
%                         ...
%                         title
%                 ...
%                 vtc_pid_settings
%
%     Each dataset typically has a different chemical composition, and therefore
%     must be processed uniquely. Further, datasets often contain different
%     experiment types (producing different data), e.g. NOESYPR1D or HSQCETGPSISP.
%     These typically require different processing schemes as well. In this
%     workflow, I separate these data early on under the sample's directory
%     structure. Each sample (s) thus contains one or more expTypes (t), where
%     sample(s).expType(t) indexes a "specList". A specList is a continuous
%     time series on a single sample with made up of only one data type. In
%     other words, a specList is the largest group of spectra that can be
%     processed as together in a meaningful way.
%
%     The primary steps that we typically take to converting these Bruker
%     files to MATLAB-readable data are:
%
%         1) Getting the data off the magnet and onto our local machine/Dropbox
%             In my hands, this is best done by compressing the datasets into
%             a zip/tar file on the magnet, then uploading to Dropbox.
%         2) Unpacking the data and making sure it's organized and TRACEABLE back
%             to the magnet's data storage/backups
%             I put the zip/tar file in the appropriate project directory, and
%             make a
%         3) Re-organizing the data into a directory structure that suits analysis
%             goals and reproducibility standards
%         4) Pre-processing the data using nmrPipe (there are several reasons
%             we use this software over others, including:
%                 - the software is open-source,maintained by NIST, and ubiquitous
%                 - the software is wildly flexible, and offers hundreds of options
%                 - quality of output is typically superior in our hands
%             	- easy to call from/interface with other pipelines
%         5) Loading the data into a MATLAB workflow pipeline
%             NOTE: at this point, data from related datasets may be combined
%         6) Processing the data within MATLAB using Edison Lab Metabolomics Toolbox
%             (https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA)
%         7) Visualizing the dataset (often in ways that would make a cool
%             T-shirt or album cover)
%         8) Extracting quantifiable features from the dataset.
%
%     This function covers steps 2-4. It first divides up the data and
%     organizes it into specLists sampele(s).expType(t) under a studyInfo
%     structure. It builds out the standard directory structure, moves the
%     data accordingly (to eliminate human errors/save TONS of time), and then
%     pre-processes the data in an interactive way.
%
%       Specifically, for each spectrum in each specList, we need to:
%           - produce the fid.com file
%           - run the fid.com file to produce an .fid file
%           - produce an ft.com file
%           - run the ft.com file to produce an .ft file
%           - Visualize at least one spectrum using nmrDraw to assess
%             p0 and p1 phasing, and potentially baseline, solvent suppression,
%             and line broadening issues
%           - Apply these corrections to all spectra in the current specList,
%             and reprocess them.
%           - Check the spectra again in nmrDraw to verify that artifacts have
%             been removed and no more have been introduced.
%
%     The next step is to pass the studyInfo to an .m file generation script
%     which will produce a basic MATLAB processing pipeline that results in
%     visualized CIVM data for each specList. Finally, the data should be
%     combined and extracted.
%
% Inputs:
%
%     studyInfo       structure obtained from running constructHRMASDirectory_4()
%
%     Optional Name-Value pairs:
%
%         varargin
%                 'repSpecName'   pass the name of the representative spectrum,
%                                 e.g. '1' or '22'
%                                 (default is first spectrum, natural number sorted)
%                 'auto' whether automatic run the program without human interaction. True or False. Default False.
% Outputs:
%
%     updatedStudyInfo    updated studyInfo structure (a few things are added,
%                         such as proc file info and file paths)
%     (processed data)    fid.com, fid, ft_com, ft files all generated and
%                         stored in their locations.
%
%
%
%     MTJ DEC2020

    
%% Handle options

    useProvidedRefSpec = 0;
    lastFTcomTemp = [];
    auto=false;
    if ~isempty(varargin)
        
        % Check for repSpecName argument
        
            % [~,ind] = ismember('repSpecName',varargin);
            ind=find(strcmp('repSpecName',varargin));
            if any(ind)
                providedRepSpecName = varargin{ind+1};
                useProvidedRefSpec = 1;
            end
            
        % Pass full path of existing template ft.com file
        
            % [~,ind] = ismember('useFT.com',varargin);
            ind=find(strcmp('useFT.com',varargin));
            if any(ind)
                knownftdir = varargin{ind+1};
            end
            
            ind=find(strcmp('auto',varargin));
            if any(ind)
                auto = varargin{ind+1};
            end
    end

    
    
%% General workflow

    for s = 1:length(studyInfo.sample)
        for t = 1:length(studyInfo.sample(s).expType)
            
            %% I want one object to access here, with everything I need

                specList = studyInfo.sample(s).expType(t);
                
            %% Process/Reprocess the data, or Skip this set?
            
                % Check to see if the ft files already exist for this
                % dataset.
                
                % If the ft files do not exist, offer to process or skip:
                                
                    if isempty(dir([specList.paths.ft,'/*.ft']))
                        
                        if auto
                          response=1;
                        else
                          response = menu(['Process ',specList.type,' spectra for sample ',studyInfo.sample(s).name,' ?'],...
                               'Process',...
                               'Skip');
                        end
                        if response ~= 1
                            continue
                        end
                        
                % If they do exist, offer to reprocess or skip:
                    else
                        
                        response = menu(['Re-process ',specList.type,' spectra for sample ',studyInfo.sample(s).name,' ? Current processing files will be copied and stored just in case.'],...
                             'Re-process',...
                             'Skip');
                        
                        % If reprocessing, copy the existing data for storage:
                        
                            if response == 1
                                % Move the files to another directory
                                    movefile(specList.paths.templates,...
                                             [specList.paths.scripts,'/',specList.type,'/previousTemplates_',   regexprep(   char(datetime(now,'ConvertFrom','datenum')), '\W','-')]);

                                % mkdir the nmrPipe_templates folder again
                                    % Make sure this is actually a dir
                                        mkdir(specList.paths.templates)
                                        mkdir([specList.paths.templates,'/representative_spectrum'])
                            else
                                continue
                            end
                         
                    end
           
            %% Go get the template processing files
                
                copyfile([templateDir,'/proc_civm.com'],specList.paths.templates )
                copyfile([templateDir,'/runFIDfiles.com'],specList.paths.templates)
                copyfile([templateDir,'/generateFTcoms.sh'],specList.paths.templates)
                copyfile([templateDir,'/runFTfiles.com'],specList.paths.templates)

            %% Check to make sure the data exist
                
            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Optimize the representative spectrum

    % We need template fid.com and ft.com files for the current dataset
                
            % Make the representative spectrum directory
            
                cd([specList.paths.templates,'/representative_spectrum'])

                % Default representative is first spectrum (natural number sort):
                
                    if ~useProvidedRefSpec
                        repSpecName = num2str(min(cellfun(@str2num,{studyInfo.sample(s).expType(t).files.name})));
                    else
                        repSpecName = providedRepSpecName;
                    end
                                        
                    specList.repSpecName = repSpecName;
                    
                specList.paths.rep_spec = [specList.paths.templates ,'/representative_spectrum/',repSpecName]; % doesn't exist yet *** FLAG ***
                mkdir(specList.paths.rep_spec);                                                                % now it does
                
            % Copy the representative spectrum to the templates dir
                
                copyfile([specList.paths.raw,'/',repSpecName],...
                            specList.paths.rep_spec)
                
            % Generate the fid.com file template for that spectrum
                cd(specList.paths.rep_spec)
                fidTemplate = getfidDotComFile('auto',auto);
                    specList.fidComTemplate = fidTemplate;
                    
             % Run the fid.com to generate test.fid
            
                !fid.com

      %% Make an fid.com file for each spectrum and run it
      
            specList.fidComFiles = makeFIDcomFiles(specList,fidTemplate);

%%
%             % Make the ft.com file from the test.fid file and template, or
%             % retrieve it from elsewhere
%
%             % Give option to pick a .com file or generate de novo using
%             basicFT1.com, or use last file
                  if exist('knownftdir','var')
                    [lastFTcom,lastFTcom_path,~,cancel] = ft_comTemplate(specList,repSpecName,lastFTcomTemp,knownftdir);
                  else
                    [lastFTcom,lastFTcom_path,~,cancel] = ft_comTemplate(specList,repSpecName,lastFTcomTemp);
                  end
                      lastFTcomTemp = [lastFTcom_path,'/',lastFTcom];
                    if cancel
                        fprintf('processCIVMdata() was cancelled prematurely')
                        updatedStudyInfo = studyInfo;
                        return
                    end
                  
      % Interface with template ft.com file (generated by basicFT1.com or inherited)
      
            cd(specList.paths.templates)
            specList.ftComTemplate = convertFTcomToTemplate(['./representative_spectrum/',repSpecName],'_ft.com',repSpecName,specList);
                    
      % Make and run the ft_com files
      
            [output] = makeAndRunFT_comFiles('generateFTcoms.sh',...   % name of generator file to be created and run
                                                'runFTfiles.com',...   % name of runner file to be created/created + run
                                                specList,'batchMode'); % don't run each file individually (run as bash loop)

            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % It is necessary at this point to keep track of what the ft file is
    % that we are operating on.
    
%         currentFTcomFile = specList.ftComTemplate.fullpath;
        

      %% Fine-tune processing on representative spectrum (while loop)
      
            doExit = false;
            doCancel = false;
            autoi=0;
            while ~or(doExit,doCancel)

                cd(specList.paths.templates)

                % Present the user with a menu of options
                if auto
                  if autoi==0
                    response=5;
                    autoi=autoi+1;
                  else
                    response=8;
                    doExit=true;
                  end
                  
                else
                        [response,doExit,doCancel] = CIVM_processingMenu();
                end

                switch response

                    % Examine the data (starting from representative spectrum)

                    case 1

                            % Option - Ask which spectrum?
                            cd(specList.paths.ft)
                            system(['nmrDraw ',repSpecName,'.ft'])

                    % Enter Phasing Corrections
                    case 2
%%
                            % Update phasing in the template file

%                                 updatePhasing_proc_civm(specList.paths.templates ,'proc_civm.com',repSpecName);
                                updatePhasing_repSpec(specList.paths.rep_spec,repSpecName); % at this point, template_ft.com should be what we work off of.

                            % Produce the ft.com file for the
                            % representative spectrum (template for the
                            % downstream steps)
                                cd(specList.paths.templates)
                                specList.ftComTemplate(end+1) = convertFTcomToTemplate(['./representative_spectrum/',repSpecName],'_ft.com',repSpecName,specList);
                                lastFTcomTemp = specList.ftComTemplate.fullpath;
                                
                            % Apply to dataset

                                output(end+1) = makeAndRunFT_comFiles('generateFTcoms.sh','runFTfiles.com',specList,'batchMode');

                    % Pass auto referencing args
                    case 3
%                             % Update phasing in the template file but only
%                             % run for ref spec
%
% %                                 updatePhasing_proc_civm(specList.paths.templates ,'proc_civm.com',repSpecName);
%                                 updatePhasing_repSpec(specList.paths.rep_spec,repSpecName); % at this point, template_ft.com should be what we work off of.
%
%                             % Produce the ft.com file for the
%                             % representative spectrum (template for the
%                             % downstream steps)
%                                 cd(specList.paths.templates)
%                                 specList.ftComTemplate(end+1) = convertFTcomToTemplate(['./representative_spectrum/',repSpecName],'_ft.com',repSpecName,specList);
%                                 lastFTcomTemp = specList.ftComTemplate.fullpath;
%
%                             % Apply to dataset
%
%                                 output(end+1) = makeAndRunFT_comFiles('generateFTcoms.sh','runFTfiles.com',specList,'batchMode',);
%

                    % Pass baseline args
                    case 4
%                         finder(lastFTcomTemp)

                    % Apply processing to dataset
                    case 5

                        output(end+1) = makeAndRunFT_comFiles('generateFTcoms.sh','runFTfiles.com',specList);

                    % Reset phasing
                    case 6

% %                             % Update phasing
% %
% %                                 updatePhasing_proc_civm(specList.paths.templates ,'proc_civm.com',repSpecName,'reset');
% %
%                             % Produce the ft.com template
%
%                                 specList.ftComTemplate(end+1) = convertFTcomToTemplate(['./representative_spectrum/',repSpecName],'_ft.com',repSpecName,specList);
%
%                             %  Reprocess the spectra
%
%                                 output(end+1) = makeAndRunFT_comFiles('generateFTcoms.sh','runFTfiles.com',specList,'batchMode');
                    
                    % Try different ft.com template
                    case 7
                        
                      % Get ftcomtemplate
                          [~,~,~,cancel] = ft_comTemplate(specList,repSpecName);
                          if cancel
                              fprintf('\n\tprocessCIVMdata() was cancelled prematurely\n\n')
                              return
                          end

                      % Interface with template ft.com file (generated by basicFT1.com or inherited)

                            cd(specList.paths.templates)
                            specList.ftComTemplate = convertFTcomToTemplate(['./representative_spectrum/',repSpecName],'_ft.com',repSpecName,specList);

                      % Make and run the ft_com files

                            [output] = makeAndRunFT_comFiles('generateFTcoms.sh',...   % name of generator file to be created and run
                                                                'runFTfiles.com',...   % name of runner file to be created/created + run
                                                                specList,'batchMode'); % don't run each file individually (run as bash loop)
                    case 8
                        fprintf('\nExiting and saving results in studyInfo struct\n\n')
                        
                      % Make sure we store all the info we've generated for this specList
            
                            specList.ftFileInfo = output;
                            studyInfo.sample(s).expType(t).paths = specList.paths;
                            studyInfo.sample(s).expType(t).type = specList.type;
                            studyInfo.sample(s).expType(t).files = specList.files;
                            studyInfo.sample(s).expType(t).repSpecName = specList.repSpecName;
                            studyInfo.sample(s).expType(t).fidComTemplate = specList.fidComTemplate;
                            studyInfo.sample(s).expType(t).ftComTemplate = specList.ftComTemplate;
                            studyInfo.sample(s).expType(t).fidComFiles = specList.fidComFiles;
                            studyInfo.sample(s).expType(t).ftFileInfo = specList.ftFileInfo;

                    case 10
                        fprintf('\nExiting without saving\n\n')
                        
                       
                end

            end
                            
        end
    end

    updatedStudyInfo = studyInfo;
%     cd(destinationDir)
end
