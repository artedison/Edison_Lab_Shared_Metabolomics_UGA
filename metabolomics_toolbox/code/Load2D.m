function spectra=Load2D(path)

%%NOTE:
% spectra = Load2D(path)
%
% Load Bruker NMR data from the processed 2rr file at 'path'.
%
% Arguments:
% path                 Path to the 2rr files to be loaded.
%
% Return Values:
% in the form of spectra.*
% real               Real part of the spectra.
% ppm1               F2 PPM scale.
% ppm2               F1 PPM scale.
%
% LOG:
% 
% Last Revision 1/5/2010
% Steven L Robinette, University of Florida
% RT applied changes from Load1D to Load2D - 28 Sept 2017
% RT change deblank to strtrim to get rid of leading and trailing
% whitespaces 30 Jan 2018
%% Load spectra
if iscell(path)==1
    label=path;
    for k=1:length(path)
        spectra(k)=Load2rr(path{k});
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
        %c=char(c);
    end
    c = sort(cellfun(@str2num,c(find(~(cellfun(@isempty,regexp(c,'^\d+$')))))));
    c = num2str(c);
    cd(a);
    m=0;
    for k=1:size(c,1)   
        if exist(strcat(path,'/',strtrim(c(k,:)),'/pdata/1/2rr'))==2
            m=m+1;
            %spectra(m)=Load2rr(strcat(path,'/',deblank(c(k,:))));
            spectra(m)=Load2rr(strcat(path,'/',strtrim(c(k,:))));
            if ~ischar(spectra(m).Title)
                spectra(m).Title='label your samples next time!';
            end
        end
    end
end

