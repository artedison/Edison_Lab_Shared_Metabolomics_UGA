function [newmFilename] = makeMfile(template,sampleName,expType,destination)
      % Edit the m file
        % Change all instances of 'sampleName' to the actual sample name
            filedata = fileread(template);
            filedata = regexprep(filedata,'sampleNameGoesHere',sampleName);
            filedata = regexprep(filedata,'expName',expType);
            filedata = regexprep(filedata,'mFileLocation',destination);           

        % Write a new, renamed file to template directory (will be moved)
            %newmFilename = regexprep(template,'sampleName',[sampleName,'_',expType]);
            cd(destination)
            newmFilename = [expType,'_short.m'];
            
            f = fopen(newmFilename,'w');
                fprintf(f,'%s',filedata);
            fclose(f);

            cd(destination)
end