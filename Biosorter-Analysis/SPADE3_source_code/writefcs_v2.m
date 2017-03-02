function writefcs(filename, data, marker_names,channel_names, fcshdr)
% writefcs(filename, data, marker_names,channel_names, fcshdr)

if exist('fcshdr') && isstruct(fcshdr)
    % write to fcs file the data, and the scaling + spillover in fcshdr
    if size(data,2)~= length(marker_names) % put the data matrix back to what flow people are familiar with, thin tall matrix
        if size(data,1)== length(marker_names)
            data = data.';
        else
            error('data size and marker_names length do not match!!')
        end
    end
    fcsheader_main=['\$BEGINANALYSIS\0\$ENDANALYSIS\0\$BEGINSTEXT\0\$ENDSTEXT\0\$NEXTDATA\0\'];
    fcsheader_main = [fcsheader_main,'$TOT\',num2str(size(data,1)),'\']; % number of cells/events
    fcsheader_main = [fcsheader_main,'$PAR\',num2str(size(data,2)),'\']; % number of channels
    fcsheader_main = [fcsheader_main,'FCSversion\3\'];  % i'm pretending this is a fcs v3 format
    fcsheader_main = [fcsheader_main,'CREATOR\PengQiu FCS writer\'];  
    [PATHSTR,NAME,EXT] = fileparts(filename); fcsheader_main = [fcsheader_main,'FILENAME\',[NAME,EXT],'\'];  
    fcsheader_main = [fcsheader_main,'$BYTEORD\4,3,2,1\'];  % don't know what this means
    fcsheader_main = [fcsheader_main,'$DATATYPE\F\'];  % float precision, 32 bits
    fcsheader_main = [fcsheader_main,'$MODE\L\'];  % list mode
    for i=1:length(marker_names)
        fcsheader_main = [fcsheader_main,'$P',num2str(i),'B\',num2str(32),'\'];
        if exist('channel_names')
            fcsheader_main = [fcsheader_main,'$P',num2str(i),'N\',channel_names{i},'\'];
        else
            fcsheader_main = [fcsheader_main,'$P',num2str(i),'N\',marker_names{i},'\'];
        end
        fcsheader_main = [fcsheader_main,'$P',num2str(i),'S\',marker_names{i},'\'];
        
        fcsheader_main = [fcsheader_main,'$P',num2str(i),'E\',num2str(fcshdr.par(i).decade),',',num2str(fcshdr.par(i).logzero),'\'];
        if isequal(fcshdr.par(i).G,1) || isequal(fcshdr.par(i).G,0) || isempty(fcshdr.par(i).G) || fcshdr.par(i).log~=0
            1; % do nothing
        else
            fcsheader_main = [fcsheader_main,'$P',num2str(i),'G\',num2str(fcshdr.par(i).G),'\'];
        end
        fcsheader_main = [fcsheader_main,'$P',num2str(i),'R\',num2str(fcshdr.par(i).range),'\'];
    end
    
    spillover = [num2str(length(fcshdr.compensated_channels)),cell2mat(strcat(',',channel_names(fcshdr.compensated_channels))')];
    spillover_matrix = fcshdr.spillover_matrix(fcshdr.compensated_channels,fcshdr.compensated_channels);
    for i=1:size(spillover_matrix,1)
        for j=1:size(spillover_matrix,2)
            spillover = [spillover,',',num2str(spillover_matrix(i,j),16)];
        end
    end
    fcsheader_main = [fcsheader_main,'$SPILLOVER\',spillover,'\'];  % list mode

else
    % write to fcs file the data,n and trivial scaling and spillover
    if size(data,2)~= length(marker_names) % put the data matrix back to what flow people are familiar with, thin tall matrix
        if size(data,1)== length(marker_names)
            data = data.';
        else
            error('data size and marker_names length do not match!!')
        end
    end
    
    fcsheader_main=['\$BEGINANALYSIS\0\$ENDANALYSIS\0\$BEGINSTEXT\0\$ENDSTEXT\0\$NEXTDATA\0\'];
    fcsheader_main = [fcsheader_main,'$TOT\',num2str(size(data,1)),'\']; % number of cells/events
    fcsheader_main = [fcsheader_main,'$PAR\',num2str(size(data,2)),'\']; % number of channels
    fcsheader_main = [fcsheader_main,'FCSversion\3\'];  % i'm pretending this is a fcs v3 format
    fcsheader_main = [fcsheader_main,'CREATOR\PengQiu FCS writer\'];  
    % fcsheader_main = [fcsheader_main,'$COM\PengQiu FCS writer\'];  % comment
    % fcsheader_main = [fcsheader_main,'GUID\1.fcs\ORIGINALGUID\1.fcs\'];  % don't know what this means
    [PATHSTR,NAME,EXT] = fileparts(filename); fcsheader_main = [fcsheader_main,'FILENAME\',[NAME,EXT],'\'];  
    fcsheader_main = [fcsheader_main,'$BYTEORD\4,3,2,1\'];  % don't know what this means
    fcsheader_main = [fcsheader_main,'$DATATYPE\F\'];  % float precision, 32 bits
    fcsheader_main = [fcsheader_main,'$MODE\L\'];  % list mode
    for i=1:length(marker_names)
        fcsheader_main = [fcsheader_main,'$P',num2str(i),'B\',num2str(32),'\'];
        if exist('channel_names')
            fcsheader_main = [fcsheader_main,'$P',num2str(i),'N\',channel_names{i},'\'];
        else
            fcsheader_main = [fcsheader_main,'$P',num2str(i),'N\',marker_names{i},'\'];
        end
        fcsheader_main = [fcsheader_main,'$P',num2str(i),'S\',marker_names{i},'\'];
        fcsheader_main = [fcsheader_main,'$P',num2str(i),'R\',num2str(ceil(max(data(:,i)))),'\'];
        fcsheader_main = [fcsheader_main,'$P',num2str(i),'E\','0,0','\'];
    end
    
    compensated_channels = setdiff(1:length(channel_names),[isInList(channel_names,'time'),isInList(channel_names,'event'),isInList(marker_names,'time'),isInList(marker_names,'event')]); % find non-time and -event channels
    spillover = [num2str(length(compensated_channels)),cell2mat(strcat(',',channel_names(compensated_channels))')];
    spillover_matrix = eye(length(compensated_channels));
    for i=1:size(spillover_matrix,1)
        for j=1:size(spillover_matrix,2)
            spillover = [spillover,',',num2str(spillover_matrix(i,j))];
        end
    end
    fcsheader_main = [fcsheader_main,'$SPILLOVER\',spillover,'\'];  % list mode
    
end


fid = fopen(filename,'w','b');
HeaderStart = 100;
HeaderStop = HeaderStart + length(fcsheader_main)+100-1;
DataStart = HeaderStop;
DataEnd = DataStart+prod(size(data))*4;
fcsheader_1stline  = sprintf('FCS3.0    %8d%8d%8d%8d%8d%8d',HeaderStart,HeaderStop,DataStart,DataEnd,0,0);
fcsheader_main = [fcsheader_main,'$BEGINDATA\',num2str(DataStart),'\']; 
fcsheader_main = [fcsheader_main,'$ENDDATA\',num2str(DataEnd),'\']; 
entire_header = [fcsheader_1stline, repmat(char(32),1,HeaderStart-length(fcsheader_1stline)),fcsheader_main];
entire_header = [entire_header, repmat(char(32),1,HeaderStop-length(entire_header))];
fwrite(fid,entire_header,'char');
fwrite(fid,data.','float32');
fclose(fid);


return