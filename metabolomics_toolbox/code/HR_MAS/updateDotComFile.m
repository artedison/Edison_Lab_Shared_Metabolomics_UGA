function updateDotComFile(source,oldFilename,destination,newFilename,pipepars)
    % Navigate to the file and read it
      cd(source)
      filedata = fileread(oldFilename);

    % Update the parameters
            % Temporarily remove second phasing command (for ditching
            % imaginaries)
                [phaseLineStart,phaseLineEnd,phaseLineText] = regexp(filedata,'| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \','start','end','match');
                phaseLineText = phaseLineText{:};
                filedata = [filedata(1:phaseLineStart-1),'phaseLineText',filedata(phaseLineEnd+1:end)];

        % Transcribe the other parameters from bruk2pipe (does not run if unless operating on template, which contains 'bruk2pipeParamsGoHere' ):
                block = regexprep(pipepars(contains({pipepars.name},'bruk2pipeParamsGoHere')).value,'\n  -','\n					-');
                    block = regexprep(block,'\','\\\'); % need this to cancel out escape characters with \ followed by \n:
                    filedata = strrep(filedata,'bruk2pipeParamsGoHere',block);
                
            for p = 2:length(pipepars) % go through the parameters we have, skipping the first one (bruk2pipe block)
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
            filedata = regexprep(filedata,'\\\','\'); % need this to replace \\\n with \\n in the text
            
    % Write the new file to new scripts directory
        cd(destination)
        f = fopen(newFilename,'w');
            fprintf(f,'%s',filedata);
        fclose(f);

    % Make executable 

        fileattrib(newFilename,'+x','a');

                    
end
                   