function [startTimes,dateTimes] = getRunTimes_NMR(dataDir)
%% getStartTimes
% Get timepoints from acqu files (Bruker) and return as 
% Hours since proleptic ISO start date. Based on Load1D; the rows in the
% outputs correspond to the natural number sorted filenames as in Load1D.
% This only works for Bruker files at the moment, and is very useful for
% time-series HR-MAS runs. 

% MJ 22DEC2017

        a=cd;

if ischar(dataDir)
    % Get to the files and get their names
        cd(dataDir);
        if ispc==1 % if this is a PC and not a MAC
            c=ls;
            c=c(3:end,:);
        else
            c=strread(ls,'%s');
        end
        % Sort the filenames by natural number sort:
            c = sort(cellfun(@str2num,c(find(~(cellfun(@isempty,regexp(c,'^\d+$')))))));
            c = num2str(c);
elseif iscell(dataDir)
    % Make the c list of filenames as cell array of strings
        %dataDir = {spectra.FileName}';
        c = dataDir;
end

    % If the filenames are provided (e.g. only a subset of files are wanted)
        
    
% Initialize Variables        
    startTimes = zeros(size(c,1),1);
    dateTimes = cell(size(c,1),1); 

% Loop through the files and extract the start time
    m=0;
    for i=1:size(c,1)   
            %if exist(strcat(dataDir,'/',strtrim(c(i,:)),'/acqu'))==2
            if iscell(c(i,:))
                %name = strcat(strtrim(c{i,:}),'/acqu'); % this file can
                %change if the acquisition parameters are changed in
                %TopSpin
                name = strcat(strtrim(c{i,:}),'/acqus'); % this one is static
            else 
                %name = strcat(strtrim(c(i,:)),'/acqu'); % this file can
                %change if the acquisition parameters are changed in
                %TopSpin
                name = strcat(strtrim(c(i,:)),'/acqus'); % this one is static
            end
            if exist(name)==2
                m=m+1;
                %paramsFile = fopen(strcat(dataDir,'/',strtrim(c(i,:)),'/acqu'));
                paramsFile = fopen(name);
                    % Get the seventh line (contains the datetime on Bruker
                    % files)
                        fgetl(paramsFile);
                        fgetl(paramsFile);
                        fgetl(paramsFile);
                        fgetl(paramsFile);
                        fgetl(paramsFile);
                        fgetl(paramsFile);
                    % Parse as hours 
                        dateLine = strsplit(fgetl(paramsFile),' ');
                        startDate = [dateLine{2},'-',dateLine{3}];
                        dateTimes{i} = datetime(startDate,'InputFormat','yyyy-MM-dd-HH:mm:ss.SSS'); %equal to the number of days since January 0, 0000 in the proleptic ISO calendar
                        startTimes(i) = datenum(dateTimes{i})*24; % now it's in hours            
                fclose(paramsFile);
            else fprintf(['\t No file "',strtrim(c(i,:)),'" found...\n'])
            end
    end
    cd(a);
end