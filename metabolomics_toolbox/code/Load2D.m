function spectra=Load2D(path)

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
% Last Revision 1/5/2010
% Steven L Robinette, University of Florida


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
        c=char(c);
    else
        c=strread(ls,'%s');
        c=char(c);
    end
    
    cd(a);
    
    for k=1:size(c,1)
        spectra(k)=Load2rr(strcat(path,'/',deblank(c(k,:))));
    end
end

