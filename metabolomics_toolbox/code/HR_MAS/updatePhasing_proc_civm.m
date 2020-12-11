function [ftPars,ftFile] = updatePhasing_proc_civm(pipeTempLoc,ftTempFilename)

    %pipeTempLoc = pipeTempLoc;
    %ftTempFilename = 'proc_civm.com';

        instructions = 'Phase Corrections entry';
        dims = [1,150];
        opts.Resize = 'on';
        opts.WindowStyle = 'modal';
        opts.Interpreter = 'none';
        corrs = inputdlg({'Enter/paste phasing corrections in the following boxes:    p0 correction:','p1 correction:'},instructions,dims,{'0','0'},opts);

        % Update 
            % Get the current ftdata (correction will be
            % reset, but the proc_civm.com file keeps the
            % actual phase values stored during multiple 
            % iterations of phasing
            
                ftFile = [pipeTempLoc,'/',ftTempFilename];
                ftData = fileread(ftFile);

            % Modify the phase values in the template file
                % Convert corrs to doubles
                    corrs = str2double(corrs); 


                % Extract out current phase vals (even
                % commented ones; watch out for that)
                ftPars = struct();

                % P0 Phase param
                    [ftPars.xP0.value,...
                     ftPars.xP0.startInd,...
                     ftPars.xP0.endInd] = regexp(ftData,['(?<=','xP0','[\s]+)','[\S]+'],'match','start','end');     % Extract current value

                    ftPars.xP0.replaceWith = num2str(   str2double(ftPars.xP0.value{1}) + corrs(1)   );             % do math + format as str
                    ftData = [ftData(1:ftPars.xP0.startInd(1)-1),ftPars.xP0.replaceWith,ftData(ftPars.xP0.endInd(1)+1:end)]; % replacement

                % P0 Phase param
                    [ftPars.xP1.value,...
                     ftPars.xP1.startInd,...
                     ftPars.xP1.endInd] = regexp(ftData,['(?<=','xP1','[\s]+)','[\S]+'],'match','start','end');     % Extract current value

                    ftPars.xP1.replaceWith = num2str(   str2double(ftPars.xP1.value{1}) + corrs(2)   );             % do math + format as str
                    ftData = [ftData(1:ftPars.xP1.startInd(1)-1),ftPars.xP1.replaceWith,ftData(ftPars.xP1.endInd(1)+1:end)]; % replacement

                    
                % Update (write) the file 
                    
                    f = fopen(ftFile,'w');
                        fprintf(f,'%s',ftData);
                    fclose(f);

                % Make executable 

                    fileattrib(ftFile,'+x','a');
end