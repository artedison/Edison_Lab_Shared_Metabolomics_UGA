function loadallft(pickFiles)
% Load all .ft files in the current directory. If any arguments are given,
% a multi-selection dialog is presented. Outputs to the workspace in a
% variable called 'spectra'
% pickFiles: the loading file directory

% Files will be loaded in numerical order if they start with a digit
% separated by an underscore.
%% this is from initial toolbox
%% YW add some documents 10/10/2018

if nargin==0
    filenames = dir('*.ft');
    filenames = {filenames.name};
    %filenames = sortFilenames(filenames); %% this doesn't work, MJ implemented 3JAN2018:
        filenums = regexp(filenames,'\d+','match');
        [~,i] = sort(cellfun(@str2num,[filenums{:}]));
        filenames = filenames(i);
    ftlist = cellfun(@(x) ([cd filesep x]),filenames,'uniformoutput',0);
elseif nargin==1
    [filenames,PathName] = uigetfile('*.ft','Select files...','MultiSelect','on');
    %filenames = sortFilenames(filenames); %% this doesn't work, MJ implemented 3JAN2018:
        filenums = regexp(filenames,'\d+','match');
        [~,i] = sort(cellfun(@str2num,[filenums{:}]));
        filenames = filenames(i);
        ftlist = cellfun(@(x) ([PathName x]),filenames,'uniformoutput',0);
end
for i=1:length(ftlist)
    spectra(i)=pipe2matlab(ftlist{i});
end
for i = 1:length(ftlist)
    spectra(i).FileName = filenames{i};
end
%% Display summary
if length(filenames)>=10
    filenames(1:10)'
else
    filenames'
end

assignin('caller', 'spectra', spectra)

end
