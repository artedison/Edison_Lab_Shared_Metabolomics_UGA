function spectra=Load1D(path,format)

%%% NOTE:
% RT and MJ altered this function with the following updates (28SEP2017):
%     - non-NMR folders in 'path' directory are ignored
%     - empty Titles in the .r files are now tolerated
%     - changed means of sorting (natural number sort)
%     - ignores empty "NMR data" folders
% MTJ altered again 31JAN2022:
%     - mix of numerical and string sample names tolerated, but if integers
%     only (typical), then natural number sorted. 
%
%test now that I have a github branch
% spectra=Load1D(path,format)
%
% Loads all Bruker NMR data from the processed 1r file at 'path'.
%
% Arguments:
% path                 Path to the directory with spectral directories.
% This
% format               Data format- 'bruker' is supported
% 
% Return Values:
% in the form of spectra.*
% real               Real part of the spectra.
% imaginary          Imaginary part of the spectra (may not be avaliable).
% ppm                The PPM scale.
% XTitles            the titles of the spectral directories 

if exist('format')==0
    format='bruker';
end
if iscell(path)==1 % if you want to pass a list of NMR files as a cell array to be read
    label=path;
    for k=1:length(path)
        if strcmp(format,'bruker')==1
            spectra(k)=LoadBruker(strcat(path{k},'/pdata/1/1r'));
        end
    end
else
    a=cd;
    cd(path);
    if ispc==1
        c=ls;
        c=c(3:end,:);
        %c=char(c);
    else
        c=strread(ls,'%s');
        %c=char(c); why would we do this??
    end
    
    notNumbers = cellfun(@isempty,regexp(c,'^\d+$'));
    if all(~notNumbers)
        c = sort( cellfun(@str2num,c(find( ~notNumbers )) ) );
        c = cellstr(num2str(c));
    else
        if ~isempty(c)
            % What if the filenames aren't numbers?
                d = dir('./**/pdata/1/1r');   % require 1r files 2 levels down for each sample
                
                % Extract out the sample names
                    ctemp = cell(length(d),1);
                    for f = 1:length(d)
                        s = strsplit(d(f).folder,'/');
                        ctemp(f) = s(end-2);
                    end
        else
            error('The provided directory does not contain a 1r file in the location ./**/pdata/1/1r')
        end
    end
                                    %     if isempty(str2num(c))==0
                                    %         c=num2str(sort(str2num(c)));
                                    %     else
                                    %         [e,d]=sort(c);  % this step actually changes the contents
                                    %     end
    cd(a);
%     OneDVector=zeros(size(c,1),1);
    m=0;
    if size(c,1) < size(c,2)
        c = c';
    end
    for k=1:size(c,1)   
        if strcmp(format,'bruker')==1
            if exist(strcat(path,'/',strtrim(c{k,:}),'/pdata/1/1r'))==2
                m=m+1;
                spectra(m)=LoadBruker(strcat(path,'/',strtrim(c{k,:}),'/pdata/1/1r')); % Does this file contain the expected /pdata/1/1r file and path
                if ~ischar(spectra(m).Title)
                    spectra(m).Title='label your samples next time!';
                end
                
%             else
%                 OneDVector(k)=1; % flag this file because it doesn't contain the /pdata/1/1r file we are trying to read
            end
        end            
    end
    m = 0;
    for k = 1:size(c,1)
        if exist(strcat(path,'/',strtrim(c{k,:}),'/pdata/1/1r'))==2
            m=m+1;
            spectra(m).FileName = strtrim(c{k,:});
        end
        
    end
%     spectra(find(OneDVector==1))=[]; % get rid of flagged files 

end

end
