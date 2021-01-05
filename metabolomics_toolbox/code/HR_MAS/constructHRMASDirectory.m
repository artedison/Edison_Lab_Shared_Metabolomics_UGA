function constructHRMASDirectory(parentDest,datadir,sampledir)

%% constructHRMASDirectory

    % Author: Yue Wu and Michael T. Judge
    % Version: 0.5
    % Tested on Matlab Version R2020a
    % Date: 2020
    %
    % Description:
    %       plot_spec plots a set of chemical shift values on a set of 1D
    %       spectra. It has added functionalities for display selection,
    %       coloring and displaying labels.
    %
    % Input:
    %       X: stack of 1D spectra
    %       ppm: chemical shift vector
    %       T: Table with the sample information.
    %           5 columns with following exact headers:
    %                Run_ID,Sample_ID,Sample_desc,Sample_grp,Yvec
    %
    % Output:
    %       A plot of the 1D spectra.
    %
    % Log:
    %       Ver 0.4: When plotting, does not use the black color anymore.
    %
    % Example run:
    %

    
  % filesys_cons(goaldir,datadir,sampledir)
  % this function construct the template folder for systematically arrange NMR data
  % It will copy the raw data to the folder structure
  % it will also copy useful script to the folder
  % list of data that exist
  % goaldir: the location to create the project folder
  % datadir: the location to fetch data; the folder is expected to have multiple samples which contains multiple time points
  % sampledir: the location that the sample script is from
  %%yue wu revise 10/10/2018
  %% MJ edits 13FEB2019:
%   - dynamic file naming (use the nmr data file name as the sample name)
%   - nmrpipe interfacing with .com file
%   - 

%% MJ edits DEC2020

    
  
  %%
  
  % copyfile(strcat(sampledir,'*.m'),'./');
  % copyfile(strcat(sampledir,'*.com'),'./');
  % copyfile(strcat(sampledir,'*.txt'),'./');
  % Locate the data files and their names
      dirs=dir(datadir);
        % Get the name of the sample
            %sampleName = 
      dirsname={dirs.name};
      dirsname=dirsname(~contains(dirsname,{'.','..','.DS_Store','.tar.gz'}));
      % What is the point of this?
          existsname=dir(strcat(parentDest,'sampleGroup'));
          existsname={existsname.name};
          existsname=existsname(~ismember(existsname,{'.','..'}));
      dirsname=setdiff(dirsname,existsname);
    % Create the new data directory (assuming one sample; gets one
    % directory), then jump into the new directory
    % the function will move new folders that does not exist there.

      cd(parentDest);
%       dirname = dirsname{:};
%       mkdir(dirname);
%       cd(['./',dirname]);
%%    
  % Copy them to the new location
      for i = 1:length(dirsname)
         matchind=regexp(dirsname{i},'^\.');
         % Check to see if this is actually a data directory; 
         % bail on cases where common non-directories are included
             if length(matchind)==1
                 continue;
             end
        
        % Make the individual experiment folder
            expname=char(dirsname(i));
            mkdir(expname)
            cd(strcat(parentDest,['/',expname]))
                    
        % Create the standard directories within that folder
            expdir=pwd();
            mkdir('data');mkdir('results');mkdir('scripts');
                cd('./data');mkdir('raw');mkdir('processed');              
                        cd('./processed');mkdir('nmrpipe');
                            cd('./nmrpipe');mkdir('fid');mkdir('ft');mkdir('proc_files');
                                cd('proc_files');mkdir('fid_com');mkdir('ft_com');
            
        % Copy the data into the appropriate directory
            % Check to see if these are different data types. If so, make two
            % directories and set a flag.
                % Get the acqus files (in the order they are read; NOT natural number sorted)
                    cd([expdir,'/data/raw']);
                    paramFiles = dir(cell2mat(strcat(datadir,dirsname(i),'/**/*acqus')));
                % Open each one and get its experiment type
                    for pfile = 1:length(paramFiles)
                        filedata = fileread([paramFiles(pfile).folder,'/',paramFiles(pfile).name]);
                        % Find the text inside the < > that immediately
                        % follows an instance of '##$EXP= ' followed by <...> 
                            % 'match' option
                                %expLine = regexp(filedata,'##\$EXP= <[\w*]+>','match');
                                % Find things that look like this '##$EXP= <...>' and return the matched text:  
                                %   '##\$EXP= <[\w*]+>','match' 
                                % Within the returned text, find and return
                                % the text inside <...>:
                                %   ['(?<=<)','\w*','(?=>)'],'match'
                            expType = regexp(   regexp(filedata,'##\$EXP= <[\w*]+>','match')    ,['(?<=<)','\w*','(?=>)'],'match');                                
                            paramFiles(pfile).experimentType = expType{:}{:};
                    end
                    clear('expType','pfile','filedata');
                    
                   % Find types of experiments from that list
%                         expTypes = unique({paramFiles.experimentType});
%                         [~,key] = ismember({paramFiles.experimentType},expTypes);
                            dataTypeKey = {'HSQCETGPSISP','13c1d';...
                                           'NOESYPR1D'   ,'1h1d' };
                        dataFiles = dir([datadir,dirsname{:}]);
                            dataFiles(ismember({dataFiles.name},{'.','..'})) = []; % remove dirs                                            
                        
%                     for t = 1:length(expTypes)
%                         % Find the file indices for this type
%                             %inds = find(key == t);
%                             inds = key == t;
%                         % Make the directory for this type
%                             mkdir(expTypes{t})
%                         % Copy over the files that are of this type
%                             % not the best way to do this, but it should
%                             % work:
%                             copyfile(cell2mat(strcat(datadir,dirsname(i),'/*')),expTypes{t}); % copy all files over                            
%                             cd(expTypes{t})
%                             eval(sprintf('rmdir %s',['''',strjoin({dataFiles(~inds).name},[''' ','''']),'''']   )) % delete those we don't want
% %                             for f = 1:length(inds)
% %                                 %[{dataFiles(inds).folder}',repmat({'/'},length(inds),[]),{dataFiles(inds).name}'];
% %                                 copyfile([dataFiles(inds(f)).folder,'/',dataFiles(inds(f)).name],expTypes{t});                        
% %                             end
%                     end
            % source: all files (*) in the data directory,
            % destination: ./ (present directory)
            copyfile(cell2mat(strcat(datadir,dirsname(i),'/*')),'./');
        % Copy over the scripts from the template directory
            cd(expdir);
            cd('./scripts');
        % Edit the m file
            % Change all instances of 'sampleName' to the actual sample name
                mFilename = dir(strcat(sampledir,'*.m'));
                    mFilename = mFilename.name;
                cd(sampledir)
                filedata = fileread(mFilename);
                filedata = regexprep(filedata,'sampleNameGoesHere',expname);
                
            % Write a new, renamed file to template directory (will be moved)
                newmFilename = regexprep(mFilename,'sampleName',expname);
                f = fopen(newmFilename,'w');
                    fprintf(f,'%s',filedata);
                fclose(f);
                
            % Move the new .m file to the new scripts directory 
                movefile(newmFilename,[expdir,'/scripts']);
                cd([expdir,'/scripts'])
                
%             % The .com file will just be copied (for now)
%                 copyfile(strcat(sampledir,'*.com'),'./');
%                 dotcomFilename = dir(strcat(sampledir,'*.com'));
%                 dotcomFilename = dotcomFilename.name;
            %copyfile(strcat(sampledir,'*.txt'),'./');

            
        %% Edit the nmrPipe file(s): 13C1d 
            dataType = dataTypeKey(4); % hardcode (temp)
            %if strcmp(dataType,'13c1d')
            if strcmp(dataType,'1h1d')
                % Run bruker on first timepoint
                    cd ../data/raw/2    % hardcode (temp)
                    ! bruker

                % Read the Parameters from fid.com file
                    filedata = fileread('fid.com');
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


                    for p = 1:length(pipepars)
                        tmp = regexp(filedata,['(?<=',pipepars(p).name,'[\s]+)','[\S]+'],'match');
                        tmp = tmp{:};
                        numstr = str2double(tmp); % returns nan if not just a number
                        if ~isnan(numstr)
                            pipepars(p).value = numstr;
                        else
                            pipepars(p).value = tmp;
                        end
                    end

                % Update the params in the batch .com file
                    %cd ../../../scripts
                    cd(sampledir)
                        dotcomFilename = dir(strcat(sampledir,'*.com'));
                        dotcomFilename = dotcomFilename.name;
                    filedata = fileread(dotcomFilename);
                    % Do the grpdly edit, and apply to phasing:
                        % Locate and extract the parameter
                            %[startInd,endInd,tmp] = regexp(filedata,['(?<=','-grpdly','[\s]+)','[\w*]+'],'start','end','match');
                            %regexprep(filedata,[' -grpdly','[\s]+','[\w*]+'],''); % in case it's there

                        % Apply to phase corrections
                            pipepars(end+1).name = '-p1';
                            %pipepars(end).value = pipepars(~cellfun(@isempty,strfind({pipepars.name},'-grpdly'))).value * -360;                    
                            pipepars(end).value = pipepars(contains({pipepars.name},'-grpdly')).value * -360;                    
                            pipepars(end+1).name = '-p0';
                            pipepars(end).value = 0;
                    % Temporarily remove second phasing command (for ditching
                    % imaginaries)
                        [phaseLineStart,phaseLineEnd,phaseLineText] = regexp(filedata,'| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \','start','end','match');
                        phaseLineText = phaseLineText{:};
                        filedata = [filedata(1:phaseLineStart-1),'phaseLineText',filedata(phaseLineEnd+1:end)];
                
                % Transcribe the other parameters from bruk2pipe:

                    for p = 1:length(pipepars) % go through the parameters we have
                        %[startInd,endInd,tmp] = regexp(filedata,['(?<=',pipepars(p).name,'[\s]+)','[\w*]+'],'start','end','match');
                        [startInd,endInd,~] = regexp(filedata,['(?<=',pipepars(p).name,'[\s]+)','[\S]+'],'start','end','match');
                        if ~isempty(startInd)  % if the parameter was found
                            % Before writing, all values must be strings
                            if ~isa(pipepars(p).value,'string') % assume it's a number
                                filedata  = [filedata(1:startInd-1),num2str(pipepars(p).value,8),filedata(endInd+1:end)];
                            else 
                                filedata  = [filedata(1:startInd-1),pipepars(p).value,filedata(endInd+1:end)];
                            end
                        end
                    end              

                % Add back the second phasing line
                    filedata = regexprep(filedata,'phaseLineText',phaseLineText);

                % Write the new file to new scripts directory
                    cd(expdir)
                    cd scripts
                    newFilename = regexprep(dotcomFilename,'\.com','_updated.com');        
                    f = fopen(newFilename,'w');
                        fprintf(f,'%s',filedata);
                    fclose(f);

                % Make executable 
                
                    fileattrib(newFilename,'+x','a');
                    
            % Run .com file with the number of spectra as argument

                eval(sprintf('! %s %i',newFilename,length(dataFiles)-1)) % -1 in case the last file was incomplete

                % Run nmrDraw in the nmrpipe/ft directory after previous step:

                    cd ../data/nmrpipe/ft_1h1d
                    ! nmrDraw

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
                    p0corr = str2double(corrs{1});
                    p1corr = str2double(corrs{2}); 

                % Update the .com file again and re-run nmrPipe:
                    % Read the file
                      cd(expdir)
                      cd scripts
                      filedata = fileread(newFilename);

                    % Temporarily remove second phasing command (for ditching imaginaries)
                        [phaseLineStart,phaseLineEnd,phaseLineText] = regexp(filedata,'| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \','start','end','match');
                        phaseLineText = phaseLineText{:};
                        filedata = [filedata(1:phaseLineStart-1),'phaseLineText',filedata(phaseLineEnd+1:end)];

                    % Update phasing based on user input:
                        % Find the phasing params (still stored from file
                        % write) and correct the values
%                             p1ind = ~cellfun(@isempty,strfind({pipepars.name},'-p0'));
%                             p2ind = ~cellfun(@isempty,strfind({pipepars.name},'-p1'));
                            p1ind = find(contains({pipepars.name},'-p0'));
                            p2ind = find(contains({pipepars.name},'-p1'));
                            pipepars(p1ind).value = pipepars(p1ind).value + p0corr;
                            pipepars(p2ind).value = pipepars(p2ind).value + p1corr;

                        for p = sort([p1ind,p2ind]) % go through the parameters we have
                            [startInd,endInd,~] = regexp(filedata,['(?<=',pipepars(p).name,'[\s]+)','[\S]+'],'start','end','match');
                            if ~isempty(startInd)  % if the parameter was found
                                % Before writing, all values must be strings
                                if ~isa(pipepars(p).value,'string') % assume it's a number
                                    filedata  = [filedata(1:startInd-1),num2str(pipepars(p).value,8),filedata(endInd+1:end)];
                                else 
                                    filedata  = [filedata(1:startInd-1),pipepars(p).value,filedata(endInd+1:end)];
                                end
                            end
                        end              

                % Add back the second phasing line
                    filedata = regexprep(filedata,'phaseLineText',phaseLineText);

                % Rewrite the new file
                    %newFilename = regexprep(dotcomFilename,'\.com','_updated.com');        
                    f = fopen(newFilename,'w');
                        fprintf(f,'%s',filedata); % this should overwrite the old data
                    fclose(f);      
                    
                % Re-run the updated .com file
                    eval(sprintf('! %s %i',newFilename,length(dataFiles)-1)) % -1 in case the last file was incomplete
                    
            end               

%         %% Edit the nmrPipe file(s): 1h1d
%             % Run bruker on first timepoint
%                 cd ../data/raw/2
%                 ! bruker
% 
%             % Read the Parameters from fid.com file
%                 filedata = fileread('fid.com');
%                 pipepars = struct();
%                 pipepars(1).name = '-xN';
%                 pipepars(2).name = '-xT';
%                 pipepars(3).name = '-xMODE';
%                 pipepars(4).name = '-xSW';
%                 pipepars(5).name = '-xOBS';
%                 pipepars(6).name = '-xCAR';
%                 pipepars(7).name = '-xLAB';
%                 pipepars(8).name = '-ndim';
%                 pipepars(9).name = '-grpdly';
%                 pipepars(10).name = '-bad';
%                 pipepars(11).name = '-decim';
%                 pipepars(12).name = '-dspfvs';
%                 pipepars(13).name = '-ws';
% 
% 
%                 for p = 1:length(pipepars)
%                     tmp = regexp(filedata,['(?<=',pipepars(p).name,'[\s]+)','[\w*]+'],'match');
%                     tmp = tmp{:};
%                     numstr = str2double(tmp); % returns nan if not just a number
%                     if ~isnan(numstr)
%                         pipepars(p).value = numstr;
%                     else
%                         pipepars(p).value = tmp;
%                     end
%                 end
% 
%                 %expType = regexp(       ,['(?<=<)','\w*','(?=>)'],'match');
%             % Modify the batch .com file with the parameters
%                 cd ../../../scripts
%                 filedata = fileread([sampledir,dotcomFilename]);
%                 % Do the grpdly edit, and apply to phasing:
%                     % Locate and extract the parameter
%                         %[startInd,endInd,tmp] = regexp(filedata,['(?<=','-grpdly','[\s]+)','[\w*]+'],'start','end','match');
%                         regexprep(filedata,[' -grpdly','[\s]+','[\w*]+'],''); % in case it's there
% 
%                     % Apply to phase corrections
%                         pipepars(end+1).name = '-p1';
%                         pipepars(end).value = pipepars(~cellfun(@isempty,strfind({pipepars.name},'-grpdly'))).value * -360;                    
%                         pipepars(end+1).name = '-p0';
%                         pipepars(end).value = 0;
%                 % Temporarily remove second phasing command (for ditching
%                 % imaginaries)
%                     [phaseLineStart,phaseLineEnd,phaseLineText] = regexp(filedata,'| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \','start','end','match');
%                     phaseLineText = phaseLineText{:};
%                     filedata = [filedata(1:phaseLineStart-1),'phaseLineText',filedata(phaseLineEnd+1:end)];
% 
% %                        
%             % Transcribe the other parameters from bruk2pipe:
% 
%                 for p = 1:length(pipepars) % go through the parameters we have
%                     [startInd,endInd,~] = regexp(filedata,['(?<=',pipepars(p).name,'[\s]+)','[\S]+'],'start','end','match');
%                     if ~isempty(startInd)  % if the parameter was found
%                         % Before writing, all values must be strings
%                             if ~isa(pipepars(p).value,'string') % assume it's a number
%                                 filedata  = [filedata(1:startInd-1),num2str(pipepars(p).value,8),filedata(endInd+1:end)];
%                             else 
%                                 filedata  = [filedata(1:startInd-1),pipepars(p).value,filedata(endInd+1:end)];
%                             end
%                     end
%                 end              
% 
%             % Add back second phasing line
%                 filedata = regexprep(filedata,'phaseLineText',phaseLineText);
% 
%             % Write the file
%                 newFilename = regexprep(dotcomFilename,'\.com','_updated.com');        
%                 f = fopen(newFilename,'w');
%                 fprintf(f,'%s',filedata);
%                 fclose(f);
% 
%         % Run .com file with the number of spectra as argument
%                 eval(sprintf('! %s %i',newFilename,dataFiles))
% 
%         % Run nmrDraw in the nmrpipe/ft directory after previous step:
%                 cd ../data/nmrpipe/
%             % Make pop-up:
%                 instructions = 'Phase Corrections entry from nmrDraw';
%                 dims = [1,550];
%                 opts.Resize = 'on';
%                 opts.WindowStyle = 'modal';
%                 opts.Interpreter = 'none';
%                 corrs = inputdlg({['On bruk2pipe menu, please hit "Read Parameters", ',...
%                     'then "Save Script" using defaults. Exit bruk2pipe. nmrDraw ',...
%                     'will pop up. Once you have your phase correction values, ',...
%                     'enter/paste them in the following boxes:    p1 correction:'],'p2 correction:'},instructions,dims,{'0','0'},opts);
%                 p0corr = str2num(corrs{1});
%                 p1corr = str2num(corrs{2}); 
                

        %% Open the m file for the user
            open(newmFilename)

      end
end
