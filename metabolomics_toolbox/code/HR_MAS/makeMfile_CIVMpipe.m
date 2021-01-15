function [newmFilename] = makeMfile_CIVMpipe(template,specList,sampleName)

      % Go to where the file will be placed
      
            cd([specList.paths.scripts,'/',specList.type])
            
      % Edit the m file
        % Change all instances of 'sampleName' to the actual sample name
            filedata = fileread(template);
            filedata = regexprep(filedata,'sampleNameGoesHere',sampleName);
            filedata = regexprep(filedata,'expName',specList.type);
            filedata = regexprep(filedata,'ftDataPath_relative',fullPath2RelativePath(pwd,specList.paths.ft));
            filedata = regexprep(filedata,'rawDataPath_relative',fullPath2RelativePath(pwd,specList.paths.raw));
            filedata = regexprep(filedata,'mFileLocation',specList.paths.scripts);           

        % Write a new, renamed file to template directory (will be moved)
            
            newmFilename = [specList.type,'_short.m'];
            
            f = fopen(newmFilename,'w');
                fprintf(f,'%s',filedata);
            fclose(f);

end