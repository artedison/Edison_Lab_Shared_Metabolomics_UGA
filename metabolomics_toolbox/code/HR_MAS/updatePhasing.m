function updatePhasing(filepath,oldFilename,newFilename,pipepars)

    % Navigate to the file and read it
      cd(filepath)
      filedata = fileread(oldFilename);

    % Update the parameters
            % Do the grpdly edit, and apply to phasing:
                % Locate and extract the parameter
                    %[startInd,endInd,tmp] = regexp(filedata,['(?<=','-grpdly','[\s]+)','[\w*]+'],'start','end','match');
                    %regexprep(filedata,[' -grpdly','[\s]+','[\w*]+'],''); % in case it's there

                % Apply to phase corrections
                    pipepars(end+1).name = '-p1';
                    pipepars(end).value = pipepars(contains({pipepars.name},'-grpdly')).value * -360;                    
                    pipepars(end+1).name = '-p0';
                    pipepars(end).value = 0;
            % Temporarily remove second phasing command (for ditching
            % imaginaries)
                [phaseLineStart,phaseLineEnd,phaseLineText] = regexp(filedata,'| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \','start','end','match');
                phaseLineText = phaseLineText{:};
                filedata = [filedata(1:phaseLineStart-1),'phaseLineText',filedata(phaseLineEnd+1:end)];

        % Transcribe the other parameters from bruk2pipe:

            for p = 1:length(pipepars) % go through the parameters we have
                %[startInd,endInd,tmp] = regexp(filedata,['(?<=',pipepars(p).name,'[\s]+)','[\w*]+'],'start','end','match');
                [startInd,endInd,~] = regexp(filedata,['(?<=',pipepars(p).name,'[\s]+)','[\S]+'],'start','end','match');                                                     
                if ~isempty(startInd)  % if the parameter was found
                    % Before writing, all values must be strings
                    if ~isa(pipepars(p).value,'string') % assume it's a number
                        filedata  = [filedata(1:startInd-1),num2str(pipepars(p).value,8),filedata(endInd+1:end)];
                    else 
                        filedata  = [filedata(1:startInd-1),pipepars(p).value,filedata(endInd+1:end)];
                    end
                end
            end              

        % Add back the second phasing line
            filedata = regexprep(filedata,'phaseLineText',phaseLineText);

    % Write the new file to new scripts directory
        cd([sampleDestinationDir,'/scripts'])
        %           - "[filename]_[expType].com"
        newDotcomFile = [dotcomFilename,'_',expTypes{t},'.com'];        
        f = fopen(newDotcomFile,'w');
            fprintf(f,'%s',filedata);
        fclose(f);

    % Make executable 

        fileattrib(newDotcomFile,'+x','a');
                    
                    
end
            

                    