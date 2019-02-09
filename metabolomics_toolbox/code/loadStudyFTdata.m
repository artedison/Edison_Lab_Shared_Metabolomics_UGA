function  loadStudyFTdata( csvFileName )
% Author: Edison Lab
% Version: 0.1
% Date: 30NOV2018
%
% Description:
%       Loads all .ft data files for a study from different paths specified by a separate 
%       CSV file. All .ft files should contain data of same dimensions and
%       sizes.
%       NOTE: this function is obsolete. loadallft is preferred.
%
% Input:
%       csvFileName:    The path/name of a .csv file with two columns:
%                         dir: paths to the .ft file
%                         filename: filenames for the .ft file
%                         For example, see:
%                         Dropbox (Edison_Lab@UGA)/Projects/education/Fall_2018/bioinformatics/NMR_info/1D_noesypr1d_path.csv
%
% Output:
%       spectra:        Structure containing fields:
%                             real : the 1d vector of intensities for a given spectrum
%                             ppm : the 1d vector of ppm values corresponding to each intensity in real
%                             Title : the name of the .ft file
% Log:
%       Ver 0.1 : RT, MTJ, YW annotated and checked
%
% Example run:
%
%       loadStudyFTdata('../NMR_info/1D_noesypr1d_path.csv');

% 
    FilePathColumn=1;
    FileNameColumn=2;

    ftFiles=read_mixed_csv(csvFileName,','); %Delimiter could be an input
    ftFilesList=ftFiles(2:end, FilePathColumn); % 2 is just for skipping the column names of csv file

    for i=1:length(ftFilesList)
        ftFilesList(i)=fullfile(ftFiles(i+1,FilePathColumn),ftFiles(i+1,FileNameColumn));
    end
%disp(ftFilesList)
%disp(ftFilesList{2})
%class(ftFilesList{2})
    for i=1:length(ftFilesList)
        spectra(i)=pipe2matlab(ftFilesList{i});
    end
    assignin('caller', 'spectra', spectra)
end

function lineArray = read_mixed_csv(fileName, delimiter)
  fid = fopen(fileName,'r');   %# Open the file
  lineArray = cell(100,1);     %# Preallocate a cell array (ideally slightly
                               %#   larger than is needed)
  lineIndex = 1;               %# Index of cell to place the next line in
  nextLine = fgetl(fid);       %# Read the first line from the file
  while ~isequal(nextLine,-1)         %# Loop while not at the end of the file
    lineArray{lineIndex} = nextLine;  %# Add the line to the cell array
    lineIndex = lineIndex+1;          %# Increment the line index
    nextLine = fgetl(fid);            %# Read the next line from the file
  end
  fclose(fid);                 %# Close the file
  lineArray = lineArray(1:lineIndex-1);  %# Remove empty cells, if needed
  for iLine = 1:lineIndex-1              %# Loop over lines
    lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
                        'Delimiter',delimiter);
    lineData = lineData{1};              %# Remove cell encapsulation
    if strcmp(lineArray{iLine}(end),delimiter)  %# Account for when the line
      lineData{end+1} = '';                     %#   ends with a delimiter
    end
    lineArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data
  end
end
