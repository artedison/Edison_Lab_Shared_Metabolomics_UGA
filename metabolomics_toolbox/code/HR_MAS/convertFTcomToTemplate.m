function [ftcomTemplate] = convertFTcomToTemplate(ft_com_location,ft_com_file,repSpecName,specList)

%% Make ft template

    % Read the fdata
    
        fdata = fileread([ft_com_location,'/',repSpecName,ft_com_file]);
        
    % Parameterize the fdata to prep for sed replacement & copy
    
        % Generate relative filepaths
    
           ftcom2fid = fullPath2RelativePath(specList.paths.ft_com,...
               specList.paths.fid,...
               'useEscapeCharacters');                       
           ftcom2ft = fullPath2RelativePath(specList.paths.ft_com,...
               specList.paths.ft,...
               'useEscapeCharacters'); 
           
        % Make the replacements
        
            fdata = regexprep(fdata,...
                                ['nmrPipe -in ',ft_com_location,'/test'],...
                                ['nmrPipe -in ',ftcom2fid,'/specNumber']);

            fdata = regexprep(fdata,...
                                ['-out ',ft_com_location,'/',repSpecName,'.ft -ov'],...
                                ['-out ',ftcom2ft,'/specNumber.ft -ov']);
                    
        % Write the file            
            template_ftcom = [specList.paths.templates,'/template_ft.com'];
            f = fopen(template_ftcom,'w');
                fprintf(f,'%s',fdata);
            fclose(f);        
            
        % Output 
            ftcomTemplate.fullpath = template_ftcom;
            ftcomTemplate.fdata = fdata;
            
end