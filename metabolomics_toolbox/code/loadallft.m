function loadallft(pickFiles)
% Load all .ft files in the current directory. If any arguments are given,
% a multi-selection dialog is presented. Outputs to the workspace in a
% variable called 'spectra'

% Files will be loaded in numerical order if they start with a digit
% separated by an underscore.

if nargin==0
    filenames = dir('*.ft');
    filenames = {filenames.name};
    filenames = sortFilenames(filenames);
    ftlist = cellfun(@(x) ([cd filesep x]),filenames,'uniformoutput',0);
elseif nargin==1
    [filenames,PathName] = uigetfile('*.ft','Select files...','MultiSelect','on');
    filenames = sortFilenames(filenames);
    ftlist = cellfun(@(x) ([PathName x]),filenames,'uniformoutput',0);
end
for i=1:length(ftlist)
    spectra(i)=pipe2matlab(ftlist{i});
end

%% Display summary
if length(filenames)>=10
    filenames(1:10)'
else
    filenames'
end

assignin('caller', 'spectra', spectra)