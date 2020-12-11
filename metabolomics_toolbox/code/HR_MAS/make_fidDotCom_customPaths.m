function [newName,newLocation,newFile] = make_fidDotCom_customPaths(templateFileData,indir,outdir,outdirFull,specName)
            
    % Reg exprep to add desired params
    
        % Set up to accept input and output parameters
            newFileData = regexprep(templateFileData,'./fid',[indir,'/',specName,'/fid']); % input file path
            newFileData = regexprep(newFileData,'./test.fid',[outdir,'/',specName,'.fid']);
            
    % Write the file, pass the location, etc.
        
        %cd(destination)
        newName = [specName,'_fid.com'];
        newFile = [outdirFull,'/',newName];
        newLocation = outdirFull;       
        
        f = fopen(newFile,'w'); % use the spec number as the filename
            fprintf(f,'%s',newFileData);
        fclose(f);

    % Make executable 

        fileattrib(newFile,'+x','a');
      
    
end
