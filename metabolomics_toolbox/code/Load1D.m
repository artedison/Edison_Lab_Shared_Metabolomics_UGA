function spectra=Load1D(path,format)

%test now that I have a github branch
% spectra=Load1D(path,format)
%
% Loads all Bruker NMR data from the processed 1r file at 'path'.
%
% Arguments:
% path                 Path to the directory with spectral directories
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

if iscell(path)==1
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
        c=char(c);
    else
        c=strread(ls,'%s');
        c=char(c);
    end
    
    if isempty(str2num(c))==0
        c=num2str(sort(str2num(c)));
    else
        c=sort(c);
    end
    cd(a);
    
    OneDVector=zeros(size(c,1),1);
    for k=1:size(c,1)
        
        if strcmp(format,'bruker')==1
            if exist(strcat(path,'/',strtrim(c(k,:)),'/pdata/1/1r'))==2
                spectra(k)=LoadBruker(strcat(path,'/',strtrim(c(k,:)),'/pdata/1/1r'));
            else
                OneDVector(k)=1;
            end
        end            
    end
    spectra(find(OneDVector==1))=[];
end

end
