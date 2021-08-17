function [filename,wavelengths,timepoints,data] = readHPLCfile(filename)

% https://www.mathworks.com/help/matlab/math/interpolating-gridded-data.html

    fdata = readmatrix(filename);
    
%     sample.filename = filename;
%     sample.wavelengths = fdata(1,2:end);
%     sample.timepoints = fdata(2:end,1);
%     sample.data = fdata(2:end,2:end);

    filename = filename;
    wavelengths = fdata(1,2:end);
    timepoints = fdata(2:end,1);
    data = fdata(2:end,2:end);
    
    if ~(size(data,1) == length(timepoints))
        error(['readHPLCfile: in file ''',filename, ''', timepoints (column 1) do not match data.'])
    end
    
    if ~(size(data,2) == length(wavelengths))       
        error(['readHPLCfile: in file ''',filename, ''', wavelengths (row 1) do not match data.'])
    end
        
end