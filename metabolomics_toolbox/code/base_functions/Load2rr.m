% spectra = Load2rr(path)
%
% Load Bruker NMR data from the processed 2rr file at 'path'.
%
% RT Log: Extracts titles
% Date: 28 Sept 2017
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


function spectra = Load2rr(path)

% Split the real file path
pathstr=[path,'/pdata/1/']; 

% From the 1r file path we can get the:
% 1i file path
realFilePath = [pathstr, '2rr'];
% procs file path
procs1Path = [pathstr, 'procs'];
procs2Path = [pathstr, 'proc2s'];
% title file path
titleFilePath = fullfile(pathstr, 'title');

% Now extract first dimension perameters from the procs file, modified from
% Jake T.M. Pearce, Imperial College London
procs1 = fopen(procs1Path,'rt');
while true
        % Get line
        line = fgetl(procs1);
    
        if (~isempty(findstr('##$OFFSET=',line)))
            ind = findstr(' ',line);
            procs.offset = str2double(line(ind+1:length(line)));
        end
        
        if (~isempty(findstr('##$XDIM=',line)))
            ind = findstr(' ',line);
            procs.xdim = str2double(line(ind+1:length(line)));
        end
        
        if (~isempty(findstr('##$SW_p=',line)))
            ind = findstr(' ',line);
            procs.sw = str2double(line(ind+1:length(line)));
        end
        
        if (~isempty(findstr('##$NC_proc=',line)))
            ind = findstr(' ',line);
            procs.NC_proc = str2double(line(ind+1:length(line)));
        end
        
        if (~isempty(findstr('##$SF=',line)))
            ind = findstr(' ',line);
            procs.sf = str2double(line(ind+1:length(line)));
        end
        if (~isempty(findstr('##$SI=',line)))
            ind = findstr(' ',line);
            procs.si = str2double(line(ind+1:length(line)));
        end
        if (~isempty(findstr('##$BYTORDP=',line)))
            ind = findstr(' ',line);
            procs.bytordp = str2double(line(ind+1:length(line)));
        end
        if (~isempty(findstr('##$NC_proc=',line)))
            ind = findstr(' ',line);
            procs.NC_proc = str2double(line(ind+1:length(line)));
        end
        
        % Break when we have no more text
        if ~ischar(line)
            break
        end
end

% Raise an error if something is missing
if (isempty(procs.offset) && isempty(procs.sw) && isempty(procs.sf) && isempty(procs.si) && isempty(procs.bytordp))
    fclose(procs1);
    error('Load2rr', 'Unable to load all parameters from: %s', procsFilePath); 
end
fclose(procs1);

procs.swp = procs.sw / procs.sf;
procs.dppm = procs.swp / (procs.si - 1);
spectra.ppm1 = flipud(procs.offset : -procs.dppm : (procs.offset-procs.swp));
%spectra.Title=path;

fid=fopen(titleFilePath);
spectra.Title=fgetl(fid);
fclose(fid);

% Now extract second dimension perameters from the procs file.
proc2s2 = fopen(procs2Path,'rt');
while true
        % Get line
        line = fgetl(proc2s2);
    
        if (~isempty(findstr('##$OFFSET=',line)))
            ind = findstr(' ',line);
            proc2s.offset = str2double(line(ind+1:length(line)));
        end
        
        if (~isempty(findstr('##$XDIM=',line)))
            ind = findstr(' ',line);
            proc2s.xdim = str2double(line(ind+1:length(line)));
        end
        
        if (~isempty(findstr('##$SW_p=',line)))
            ind = findstr(' ',line);
            proc2s.sw = str2double(line(ind+1:length(line)));
        end
        
        if (~isempty(findstr('##$NC_proc=',line)))
            ind = findstr(' ',line);
            proc2s.NC_proc = str2double(line(ind+1:length(line)));
        end
        
        if (~isempty(findstr('##$SF=',line)))
            ind = findstr(' ',line);
            proc2s.sf = str2double(line(ind+1:length(line)));
        end
        if (~isempty(findstr('##$SI=',line)))
            ind = findstr(' ',line);
            proc2s.si = str2double(line(ind+1:length(line)));
        end
        if (~isempty(findstr('##$BYTORDP=',line)))
            ind = findstr(' ',line);
            proc2s.bytordp = str2double(line(ind+1:length(line)));
        end
        if (~isempty(findstr('##$NC_proc=',line)))
            ind = findstr(' ',line);
            proc2s.NC_proc = str2double(line(ind+1:length(line)));
        end
        
        % Break when we have no more text
        if ~ischar(line)
            break
        end
end

proc2s.swp = proc2s.sw / proc2s.sf;
proc2s.dppm = proc2s.swp / (proc2s.si - 1);
spectra.ppm2 = flipud(proc2s.offset : -proc2s.dppm : (proc2s.offset-proc2s.swp));


% Now we have parameters load the data.
if (procs.bytordp==0) 
    machine_format = 'l'; 
else
    machine_format = 'b'; 
end

if (exist(realFilePath, 'file') == 2)
    fid = fopen(realFilePath, 'r', machine_format);

    spectrum = (fread(fid, inf, 'int32') * realpow(2, procs.NC_proc))';
    fclose(fid);
end


%Borrowed from Nils T. Nyberg, University of Copenhagen
SI1 = procs.si; SI2 = proc2s.si;
		XDIM1 = procs.xdim; XDIM2 = proc2s.xdim;

		NoSM = SI1*SI2/(XDIM1*XDIM2);    % Total number of Submatrixes
		NoSM2 = SI2/XDIM2;		 			% No of SM along F1

		spectra.real = reshape(...
				permute(...
					reshape(...
						permute(...
							reshape(spectrum,XDIM1,XDIM2,NoSM),...
						[2 1 3]),...
					XDIM2,SI1,NoSM2),...
				[2 1 3]),...
  			    SI1,SI2)';

end