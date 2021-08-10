function [ftPars,ftFile] = updatePhasing_repSpec(pipeTempLoc,ftTempFilename,repSpecName,varargin)

    %pipeTempLoc = pipeTempLoc;
    %ftTempFilename = 'proc_civm.com';
    
if contains(varargin,'reset','IgnoreCase',true)
    reset = true;
else
    reset = false;
    
    % Get the phase corrections
        instructions = 'Phase Corrections entry';
        dims = [1,150];
        opts.Resize = 'on';
        opts.WindowStyle = 'modal';
        opts.Interpreter = 'none';
        corrs = inputdlg({'Enter/paste phasing corrections in the following boxes:    p0 correction:','p1 correction:'},instructions,dims,{'0','0'},opts);
            if isempty(corrs)
                return
            end    
        % Convert corrs to doubles
            corrs = str2double(corrs);      
end


        % Update 
            % Get the current ftdata (correction will be
            % reset, but the proc_civm.com file keeps the
            % actual phase values stored during multiple 
            % iterations of phasing
            
                ftFile = [pipeTempLoc,'/',ftTempFilename];
                ftData = fileread(ftFile);


                % Extract out current phase vals (even
                % commented ones; watch out for that)
                ftPars = struct();

                % P0 Phase param
%                     [ftPars.xP0.value,...
%                      ftPars.xP0.startInd,...
%                      ftPars.xP0.endInd] = regexp(ftData,['(?<=','xP0','[\s]+)','[\S]+'],'match','start','end');     % Extract current value
                    [ftPars.xP0.value,...
                     ftPars.xP0.startInd,...
                     ftPars.xP0.endInd] = regexp(ftData,['(?<!#[^\n]*)',...     % not preceded by a # followed by 0 or more (*) non-newline chars [^\n] (e.g. not commented out)
                                            '(?<=[|] nmrPipe -fn PS -p0 )',...  % preceded by '| nmrPipe -fn PS -p0 '
                                            '\S*',...                           % expr is any number of consecutive non-whitespace chars
                                            '(?=\s)'],...                       % followed by ' '
                                            'match','start','end');             % output args
                                                         
                    if reset
                        ftPars.xP0.replaceWith = num2str(0);                                                            % do math + format as str
                    else    
                        ftPars.xP0.replaceWith = num2str(   str2double(ftPars.xP0.value{1}) + corrs(1)   );             % do math + format as str
                    end
                        ftData = [ftData(1:ftPars.xP0.startInd(1)-1),ftPars.xP0.replaceWith,ftData(ftPars.xP0.endInd(1)+1:end)]; % replacement

                % P1 Phase param
%                     [ftPars.xP1.value,...
%                      ftPars.xP1.startInd,...
%                      ftPars.xP1.endInd] = regexp(ftData,['(?<=','xP1','[\s]+)','[\S]+'],'match','start','end');     % Extract current value
                    [ftPars.xP1.value,...
                     ftPars.xP1.startInd,...
                     ftPars.xP1.endInd] = regexp(ftData,['(?<!#[^\n]*)',...     % not preceded by a # followed by 0 or more (*) non-newline chars [^\n] (e.g. not commented out)
                                            '(?<=[|] nmrPipe -fn PS -p0 \S* -p1 )',...  % preceded by '| nmrPipe -fn PS -p0 '
                                            '\S*',...                           % expr is any number of consecutive non-whitespace chars
                                            '(?=[\s]+-di)'],...                       % followed by one or more (+) spaces, then '-di'
                                            'match','start','end');             % output args
                 
                    if reset
                        ftPars.xP1.replaceWith = num2str(0);                                                            % do math + format as str
                    else    
                        ftPars.xP1.replaceWith = num2str(   str2double(ftPars.xP1.value{1}) + corrs(2)   );             % do math + format as str
                    end
                        ftData = [ftData(1:ftPars.xP1.startInd(1)-1),ftPars.xP1.replaceWith,ftData(ftPars.xP1.endInd(1)+1:end)]; % replacement

                    
                % Update (write) the file 
                    
                    f = fopen(ftFile,'w');
                        fprintf(f,'%s',ftData);
                    fclose(f);

                % Make executable 

                    fileattrib(ftFile,'+x','a');
                    
 %% Run the file (generate ft.com) and convert to template
       % From the proc_civm.com template file:
    locs = [' ./representative_spectrum/',repSpecName,'/test.fid',... % data IN location
            ' ./representative_spectrum/',repSpecName,... % ft.com directory
            ' ./representative_spectrum/',repSpecName,... % data OUT location
            ' ',repSpecName];                      % spec name (used for naming files)

%      % Re-run the proc_civm.com file
% 
%         system(['proc_civm.com',locs])    
 
 %% Re-run the the rep spectrum file

    locs
 
end