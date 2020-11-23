function [output] = constructHRMASDirectory_2(destinationDir,newDataDir,templateDir)
%%
% Note: replaced by constructHRMASDirectory.m
%
% destinationDir = goaldir;
% newDataDir = datadir;
% templateDir = sampledir;
% destinationDir = goaldir;
% clear goaldir newDataDir sampledir

% Next: don't overwrite existing files

% Define metadir
      %metadir = destinationDir;
% Get Sample Names and Locations
      newDataDirs=dir(newDataDir);
      newDataDirs = newDataDirs(~contains({newDataDirs.name},{'.','..','.DS_Store','.tar.gz'}));

% Make a Directory to hold the Sample Directories (if it doesn't exist)
      %mkdir(metaDir);
      %cd(['./',metaDir]);

% For each sample:
for s = 1:length(newDataDirs)
    % Get the Sample Name
        sampleName = newDataDirs(s).name;       
    % Make the standard new directories and name them
            cd(destinationDir)
                mkdir(sampleName)
            sampleDestinationDir = [destinationDir,'/',sampleName];
            cd(sampleDestinationDir)
                mkdir('data');
                mkdir('results');
                mkdir('scripts');
                cd('./data');
                mkdir('nmrpipe');mkdir('raw');
            
            sampleSourceDir = [newDataDir,sampleName];
           
    % Sort the Experiments into Separate Directories by Experiment Type
        
        % Generate map: <experiment types> - filenames
        
            % Get the acqus files (in the order they are read; NOT natural number sorted)
                
                paramFiles = dir([sampleSourceDir,'/**/*acqus']); % use /**/* for recursive search
                
            % Open each one and get its experiment type
            
                for pfile = 1:length(paramFiles)
                    % Read the file as a singular string
                    
                        filedata = fileread([paramFiles(pfile).folder,'/',paramFiles(pfile).name]);
                        
                    % Find the text inside the < > that immediatelyfollows an instance of '##$EXP= ' followed by <...> 
                            %expLine = regexp(filedata,'##\$EXP= <[\w*]+>','match');
                            % Find things that look like this '##$EXP= <...>' and return the matched text:  
                            %   '##\$EXP= <[\w*]+>','match' 
                            % Within the returned text, find and return
                            % the text inside <...>:
                            %   ['(?<=<)','\w*','(?=>)'],'match'
                            
                        expType = regexp(   ...
                                            regexp(filedata,'##\$EXP= <[\w*]+>','match'),...
                                            ['(?<=<)','\w*','(?=>)'],'match'); 
                    % Store in new field in the paramFiles struct                    
                        paramFiles(pfile).experimentType = expType{:}{:}; % get the type as character vector
                end
                clear('expType','pfile','filedata'); % clean up temp vars
      
        % Define the <experiment types> 
               % Find types of experiments from that list
                    expTypes = unique({paramFiles.experimentType});
                    [~,key] = ismember({paramFiles.experimentType},expTypes);
                    dataFiles = dir(sampleSourceDir);
                        dataFiles(ismember({dataFiles.name},{'.','..'})) = []; % remove unnecessary dirs                             

        % Make a directory for each <experiment type> within "raw"
            cd([sampleDestinationDir,'/data/raw'])
                for t = 1:length(expTypes)
                    % Make the directory for this type
                        mkdir(expTypes{t})
                end        
                
        % Copy the experiments for each <experiment type>
                for t = 1:length(expTypes)
                    % Find the file indices for this type
                        %inds = key == find(strcmp(expTypes,dataTypeKey));
                        inds = key == t;
                    % Use system commands to copy the files of that type
                        cd(sampleSourceDir);
                        % use escape characters for dropbox path
                            cmdLineExpdir = regexprep(sampleDestinationDir,' ','\\ ');
                            cmdLineExpdir = regexprep(cmdLineExpdir,'(','\\(');
                            cmdLineExpdir = [regexprep(cmdLineExpdir,')','\\)'),'/data/raw/',expTypes{t}];
                        %cmd = ['cp -R ', sprintf('%s ',dataFiles(inds).name), cmdLineExpdir];
                        cmd = ['mv ', sprintf('%s ',dataFiles(inds).name), cmdLineExpdir];
                        system(cmd);                    
                end          
                cd([sampleDestinationDir,'/data/raw/'])
                
        % Make fid and ft directories for each <experiment type> within "nmrpipe"
                for t = 1:length(expTypes)
                    cd([sampleDestinationDir,'/data/nmrpipe/'])
                        mkdir(expTypes{t})
                        cd(expTypes{t})
                            mkdir('fid')
                            mkdir('ft')           
                end  
                cd([sampleDestinationDir,'/data/nmrpipe/'])
                
    % Set up the initial .m and .com files

        % for each <experiment type>
            for t = 1:length(expTypes)
            % go to the first file on the list for that type
                inds = find(key == t);
                cd([sampleDestinationDir,'/data/raw/',expTypes{t},'/',dataFiles(inds(1)).name])
                %cd(paramFiles(inds(1)).folder)
            % run bruk2pipe 
                !bruker
            % check for procParm.txt (new nmrPipe output)
                templist = dir(pwd);
                if any(contains({templist.name},'procParm.txt'))
                    % Read Params
                end
            % if not, extract the necessary parameters from fid.com and store in a params sub struct		
                if any(contains({templist.name},'fid.com'))
                    % Set up to read standard parameters:                  
                        pipepars = struct();
                            pipepars(1).name = '-xN';
                            pipepars(2).name = '-xT';
                            pipepars(3).name = '-xMODE';
                            pipepars(4).name = '-xSW';
                            pipepars(5).name = '-xOBS';
                            pipepars(6).name = '-xCAR';
                            pipepars(7).name = '-xLAB';
                            pipepars(8).name = '-ndim';
                            pipepars(9).name = '-grpdly';
                            pipepars(10).name = '-bad';
                            pipepars(11).name = '-decim';
                            pipepars(12).name = '-dspfvs';
                            pipepars(13).name = '-ws';
                        
                    % Read the parameters from fid.com:
                        filedata = fileread('fid.com');

                        for p = 1:length(pipepars)
                            % Locate the parameter and extract (as a string) the value that follows
                                tmp = regexp(filedata,['(?<=',pipepars(p).name,'[\s]+)','[\S]+'],'match');
                                numstr = str2double(tmp{:}); % returns nan if not just a number
                            % Convert strings to numbers where appropriate
                            % (make sure precision is maxed out)
                                if ~isnan(numstr)
                                    pipepars(p).value = numstr;
                                else
                                    pipepars(p).value = tmp{:};
                                end
                        end                    
                else
                    fprintf('\n\n\tNo fid.com file was found. \n\n')
                end
            % Copy generic processing file, update it, write it to scripts directory 

                % Read the file
                    cd(templateDir)
                        dotcomFilename = dir('*.com');
                        dotcomFilename = dotcomFilename.name; 
                        newDotcomFile = [dotcomFilename,'_',expTypes{t},'.com'];
                        
             % Do the grpdly edit, and apply to phasing:
                % Locate and extract the parameter
                    %[startInd,endInd,tmp] = regexp(filedata,['(?<=','-grpdly','[\s]+)','[\w*]+'],'start','end','match');
                    %regexprep(filedata,[' -grpdly','[\s]+','[\w*]+'],''); % in case it's there

                % Apply to phase corrections
                    pipepars(end+1).name = '-p1';
                    pipepars(end).value = pipepars(contains({pipepars.name},'-grpdly')).value * -360;                    
                    pipepars(end+1).name = '-p0';
                    pipepars(end).value = 0;         
                                
                % Update the parameters
                    updateDotComFile(templateDir,dotcomFilename,[sampleDestinationDir,'/scripts'],newDotcomFile,pipepars)
%{
%                         % Do the grpdly edit, and apply to phasing:
%                             % Locate and extract the parameter
%                                 %[startInd,endInd,tmp] = regexp(filedata,['(?<=','-grpdly','[\s]+)','[\w*]+'],'start','end','match');
%                                 %regexprep(filedata,[' -grpdly','[\s]+','[\w*]+'],''); % in case it's there
% 
%                             % Apply to phase corrections
%                                 pipepars(end+1).name = '-p1';
%                                 %pipepars(end).value = pipepars(~cellfun(@isempty,strfind({pipepars.name},'-grpdly'))).value * -360;                    
%                                 pipepars(end).value = pipepars(contains({pipepars.name},'-grpdly')).value * -360;                    
%                                 pipepars(end+1).name = '-p0';
%                                 pipepars(end).value = 0;
%                         % Temporarily remove second phasing command (for ditching
%                         % imaginaries)
%                             [phaseLineStart,phaseLineEnd,phaseLineText] = regexp(filedata,'| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \','start','end','match');
%                             phaseLineText = phaseLineText{:};
%                             filedata = [filedata(1:phaseLineStart-1),'phaseLineText',filedata(phaseLineEnd+1:end)];
% 
%                     % Transcribe the other parameters from bruk2pipe:
% 
%                         for p = 1:length(pipepars) % go through the parameters we have
%                             %[startInd,endInd,tmp] = regexp(filedata,['(?<=',pipepars(p).name,'[\s]+)','[\w*]+'],'start','end','match');
%                             [startInd,endInd,~] = regexp(filedata,['(?<=',pipepars(p).name,'[\s]+)','[\S]+'],'start','end','match');                                                     
%                             if ~isempty(startInd)  % if the parameter was found
%                                 % Before writing, all values must be strings
%                                 if ~isa(pipepars(p).value,'string') % assume it's a number
%                                     filedata  = [filedata(1:startInd-1),num2str(pipepars(p).value,8),filedata(endInd+1:end)];
%                                 else 
%                                     filedata  = [filedata(1:startInd-1),pipepars(p).value,filedata(endInd+1:end)];
%                                 end
%                             end
%                         end              
% 
%                     % Add back the second phasing line
%                         filedata = regexprep(filedata,'phaseLineText',phaseLineText);
% 
%                 % Write the new file to new scripts directory
%                     cd([sampleDestinationDir,'/scripts'])
%                     %           - "[filename]_[expType].com"
%                             
%                     f = fopen(newDotcomFile,'w');
%                         fprintf(f,'%s',filedata);
%                     fclose(f);
% 
%                 % Make executable 
%                 
%                     fileattrib(newDotcomFile,'+x','a');
%}                        
                % Run the file
                    % Do all the navigation from here, passed as arguments:
                        cd([sampleDestinationDir,'/scripts'])
                        pathToRaw = ['../data/raw/',expTypes{t}];
                        pathToFID = ['../../nmrpipe/',expTypes{t},'/fid'];
                        pathToFT = ['../../nmrpipe/',expTypes{t},'/ft'];
                    % Build the command as a string and run it
                        cmd = [newDotcomFile,' "',pathToRaw,'" "',pathToFID,'" "',pathToFT,'"'];
                        system(cmd)
                        
                % Run NMRDraw
                    cd(['../data/nmrpipe/',expTypes{t},'/ft'])
                    !nmrDraw
                    
                % Make pop-up for user input in Matlab:

                    instructions = 'Phase Corrections entry from nmrDraw';
                    dims = [1,150];
                    opts.Resize = 'on';
                    opts.WindowStyle = 'modal';
                    opts.Interpreter = 'none';
                    corrs = inputdlg({['On bruk2pipe menu, please hit "Read Parameters", ',...
                        'then "Save Script" using defaults. Exit bruk2pipe. nmrDraw ',...
                        'will pop up. Once you have your phase correction values, ',...
                        'enter/paste them in the following boxes:    p0 correction:'],'p1 correction:'},instructions,dims,{'0','0'},opts);
                    
                % Update phasing params in pipepars struct
                    p1ind = find(contains({pipepars.name},'-p0'));
                    p2ind = find(contains({pipepars.name},'-p1'));
                    pipepars(p1ind).value = pipepars(p1ind).value + str2double(corrs{1});
                    pipepars(p2ind).value = pipepars(p2ind).value + str2double(corrs{2}); 
                        
                    % store this info in a dotcomFiles struct
                        % filename
                        % location
                        
                % Update the .com file again using pipepars
                    updateDotComFile([sampleDestinationDir,'/scripts'],newDotcomFile,[sampleDestinationDir,'/scripts'],newDotcomFile,pipepars);
                    
                % Re-run .com file:   
                    % Do all the navigation from here, passed as arguments:
                        cd([sampleDestinationDir,'/scripts'])
                        pathToRaw = ['../data/raw/',expTypes{t}];
                        pathToFID = ['../../nmrpipe/',expTypes{t},'/fid'];
                        pathToFT = ['../../nmrpipe/',expTypes{t},'/ft'];
                    % Build the command as a string and run it
                        cmd = [newDotcomFile,' "',pathToRaw,'" "',pathToFID,'" "',pathToFT,'"'];
                        system(cmd)    
                  % Run NMRDraw
%                     cd(['../data/nmrpipe/',expTypes{t},'/ft'])
%                     !nmrDraw              
        % make an m file for the experiment type
            mFilename = dir([templateDir,'/*.m']);
                mFilename = [mFilename.folder,'/',mFilename.name];
                newmFilename = makeMfile(mFilename,sampleName,expTypes{t},[sampleDestinationDir,'/scripts']);
            cd(sampleDestinationDir),cd scripts
        % store this in an mFiles struct
            % filename
            % location
        % Open up the .m file
            cd(sampleDestinationDir),cd scripts
            open(newmFilename)

            end
    % Run the .com files
        % For each 
            
    % For a Given Data Type (Presented in Menu)
    
    % Look at the Data in nmrPipe, do Phasing
    
    % Re-run the data in nmrPipe
    
end
    output = struct();

end