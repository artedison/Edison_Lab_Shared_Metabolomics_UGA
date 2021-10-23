function paramFiles = readAcqusFiles(dataDir)
%% readAcqusFiles

    % Author: Michael T. Judge
    % Version: 0.1
    % Tested on Matlab Version R2020a
    % Date: 2DEC2020
    %
    % Description:
    %       Extracts data from acqus files in a CIVM run (or any other data
    %       set). Returned as a structure. Mainly used to extract the
    %       EXPTYPE from each file, but other params can be added in the
    %       future as needed.
    %
    % Input:
    %       dataDir:    full path to directory containing the data
    %
    % Output:
    %       paramFiles: structure containing the full text of the acqus
    %                   file for each measurement (e.g. 1:n), as well as a
    %                   few key parameters in their respective formats.
    %
    % Log:
    %       Added as part of CIVM dataset pipeline rewrite. 
    %
    
%%


    paramFiles = dir([dataDir,'/**/*acqus']); % use /**/* for recursive search
                
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
                            
%                         expType = regexp(   ...
%                                             regexp(filedata,'##\$EXP= <[\w*]+>','match'),...
%                                             ['(?<=<)','\w*','(?=>)'],'match'); 
%                         expType = regexp(   ...
%                                             regexp(filedata,'##\$PULPROG= <[\w*]+>','match'),...
%                                             ['(?<=<)','\w*','(?=>)'],'match'); 
                        % The above regexp is too strict and didn't allow
                        % things like 'hsqcetgpsisp2.2' because of the '.'.
                        % The following takes any characters except '<' and
                        % '>' (one or more instances). 
                        expType = regexp(   ...
                                            regexp(filedata,'##\$PULPROG= <[^<>]+>','match'),...
                                            ['(?<=<)','[^<>]+','(?=>)'],'match'); 
                        % As a debugging note:
                            if isempty(expType)
                                warning(['No PULPROG was found in acqus file ',[paramFiles(pfile).folder,'/',paramFiles(pfile).name]])
                            end
                                        
                    % Store in new field in the paramFiles struct                    
                        paramFiles(pfile).experimentType = expType{:}{:}; % get the type as character vector
                end
                
end