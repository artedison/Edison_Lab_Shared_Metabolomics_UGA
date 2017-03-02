function spectra = LoadBruker(path, varargin)

% Split the real file path
[pathstr, name, ext] = fileparts(path); %#ok<NASGU>

% From the 1r file path we can get the:
% 1i file path
imaginaryFilePath = fullfile(pathstr, '1i');
% procs file path
procsFilePath = fullfile(pathstr, 'procs');
% title file path
titleFilePath = fullfile(pathstr, 'title');

% Now extract perameters from the procs file.
procs = fopen(procsFilePath,'rt');
while true
    % Get line
    line = fgetl(procs);
    
    if (~isempty(findstr('##$OFFSET=',line)))
        ind = findstr(' ',line);
        offset = str2double(line(ind+1:length(line)));
    end
    
    if (~isempty(findstr('##$SW_p=',line)))
        ind = findstr(' ',line);
        sw = str2double(line(ind+1:length(line)));
    end
    
    if (~isempty(findstr('##$NC_proc=',line)))
        ind = findstr(' ',line);
        NC_proc = str2double(line(ind+1:length(line)));
    end
    
    if (~isempty(findstr('##$SF=',line)))
        ind = findstr(' ',line);
        sf = str2double(line(ind+1:length(line)));
    end
    if (~isempty(findstr('##$SI=',line)))
        ind = findstr(' ',line);
        si = str2double(line(ind+1:length(line)));
    end
    if (~isempty(findstr('##$BYTORDP=',line)))
        ind = findstr(' ',line);
        bytordp = str2double(line(ind+1:length(line)));
    end
    if (~isempty(findstr('##$NC_proc=',line)))
        ind = findstr(' ',line);
        NC_proc = str2double(line(ind+1:length(line)));
    end
    
    % Break when we have no more text
    if ~ischar(line)
        break
    end
end
% Raise an error if something is missing
if (isempty(offset) && isempty(sw) && isempty(sf) && isempty(si) && isempty(bytordp))
    fclose(procs);
    error('LoadBruker:', 'Unable to load all parameters from: %s', procsFilePath);
end
fclose(procs);

% Now we have parameters load the data.
if (bytordp==0)
    machine_format = 'l';
else
    machine_format = 'b';
end

if (exist(path, 'file') == 2)
    fid = fopen(path, 'r', machine_format);
    
    spectra.real = (fread(fid, inf, 'int32') * realpow(2, NC_proc))';
    fclose(fid);
end

if (exist(imaginaryFilePath, 'file') == 2)
    fid = fopen(imaginaryFilePath, 'r', machine_format);
    
    spectra.imaginary = (fread(fid, inf, 'int32') * realpow(2, NC_proc))';
    fclose(fid);
end

% and caluculate the ppm scale.
swp = sw / sf;
dppm = swp / (si - 1);
spectra.ppm = flipud(offset : -dppm : (offset-swp));

fid=fopen(titleFilePath);
spectra.Title=fgetl(fid);
fclose(fid);

end