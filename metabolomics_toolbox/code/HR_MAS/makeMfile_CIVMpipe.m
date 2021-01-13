function [newmFilename] = makeMfile_CIVMpipe(template,specList)
    template = [templateDir,'/sampleName_short_CIVMpipe.m'];
    specList = studyInfo_rpom_1.sample(1).expType(1);
    sampleName = studyInfo_rpom_1.sample(1).name;
    
      % Edit the m file
        % Change all instances of 'sampleName' to the actual sample name
            filedata = fileread(template);
            filedata = regexprep(filedata,'sampleNameGoesHere',sampleName);
            filedata = regexprep(filedata,'expName',specList.type);
            filedata = regexprep(filedata,'mFileLocation',specList.paths.scripts);           

        % Write a new, renamed file to template directory (will be moved)
            %newmFilename = regexprep(template,'sampleName',[sampleName,'_',expType]);
            cd(destination)
            newmFilename = [expType,'_short.m'];
            
            f = fopen(newmFilename,'w');
                fprintf(f,'%s',filedata);
            fclose(f);

            cd(destination)
end