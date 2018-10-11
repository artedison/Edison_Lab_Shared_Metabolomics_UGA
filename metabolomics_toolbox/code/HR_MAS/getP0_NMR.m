function [P0] = getP0_NMR(dataDir)
%% getP0_NMR

% Author: Michael T. Judge
% Version: 0.1
% Date: 04/12/2018

% Description:
%         Get pulsewidth for noesypr1d experiments from the static acqus
%         file. This only works for Bruker files at the moment. Very useful 
%         for time-series HR-MAS runs. Modified from Load1D.m
%
% Input: 
%       dataDir: the location of the NMR data files as an overarching
%       directory that contains the files '1','2','3',...,'n'.
%
% Output: 
%       The single value of the pulse width for the spectra in µs. 
%
% Log:
%
% Example run: 
%       dataDir = data/NMR/HRMAS_ncrassa_paper_Sample_4;
%       [P0] = getP0_NMR(dataDir)

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
        c = dataDir;
end

% Go to one file and extract the P0
    m=0;
    for i=1   
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
                paramsFile = fopen(name);
                    c=textread(name,'%s','delimiter','\n');
                    Pline = strsplit(c{find(~cellfun(@isempty,strfind(c,'##$P=')))+1},' ');
                    P0 = str2num(Pline{1}); % this is the pulse width in µseconds
                fclose(paramsFile);
            else fprintf(['\t No file "',strtrim(c(i,:)),'" found...\n'])
            end
    end
    cd(a);
end