function importXCMS(filename)

% importXCMS(filename)
%
% Import .xls file of XCMS processed MS data.  
%
% Arguments:
% 
% filename        Path to .xls file containing XCMS output data
%
% Outputs:
% msmatrix        N x P matrix of raw MS data - abundances of each feature
%                 P for all N spectra processed by XCMS
% massvect        Median m/z value for each feature P
% RTvect          Median retention time for each feature P
% metadata        XCMS calculated metadata, {fold change, t-stat, p-value,
%                 median m/z, min m/z, max m/z, median RT, min RT, max RT, number of
%                 spectra with feature, number of class 1 spectra with feature, number of class 2
%                 spectra with feature

 
[num,txt,raw] = xlsread(filename);
 
 msmatrix=cell2mat(raw(2:end,15:size(raw,2)))';
 RTvect=cell2mat(raw(2:end,8))';
 massvect=cell2mat(raw(2:end,5))';
 metadata=cell2mat(raw(2:end,2:13))';
 
assignin('base', 'msmatrix', msmatrix);
assignin('base', 'RTvect', RTvect);
assignin('base', 'massvect', massvect);
assignin('base', 'metadata', metadata);