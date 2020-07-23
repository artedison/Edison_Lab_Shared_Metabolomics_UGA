function [spectra] = mnova2matlab(dataDir)

%%  mnova2matlab
%
%     Loads NMR spectra from txt files exported in MNova. Support for PCs and Macs included.
%     
%     
%     Inputs: 
%             dataDir     directory where the txt files (one for each spectrum) are stored.
%                         No other txt files should be in there. 
%     
%     Outputs: 
%             spectra     structure analogous to that produced by loadallft() and Load1D()
%     
%
%     MTJ 2020
    
%%
    cd(dataDir)
    d = dir('*.txt');
    filenames = {d.name}';
    if ~isempty(filenames)
        spectra = struct();
        if ispc
            for i = 1:length(filenames)
                filedat = readmatrix(filenames{i},'Delimiter',{'\t'},'LineEnding','\r\n');
                if filedat(1,1) < filedat(end,1)
                    spectra(i).real = flip(filedat(:,2))';
                    spectra(i).ppm = flip(filedat(:,1))';
                else
                    spectra(i).real = filedat(:,2)';
                    spectra(i).ppm = filedat(:,1)';
                end
                spectra(i).imaginary = nan;
                spectra(i).Title = filenames{i};
                spectra(i).FileName = filenames{i};
            end
        else

            for i = 1:length(filenames)
                filedat = readmatrix(filenames{i},'Delimiter',{'\t'},'LineEnding','\n');
                if filedat(1,1) < filedat(end,1)
                    spectra(i).real = flip(filedat(:,2))';
                    spectra(i).ppm = flip(filedat(:,1))';
                else
                    spectra(i).real = filedat(:,2)';
                    spectra(i).ppm = filedat(:,1)';
                end
                spectra(i).imaginary = nan;
                spectra(i).Title = filenames{i};
                spectra(i).FileName = filenames{i};
            end
        end
    else
        error('mnova2matlab: No .txt files were found in this directory')
    end
end