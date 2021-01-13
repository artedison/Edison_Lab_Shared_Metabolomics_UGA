function [file,path,cancel,success] = ft_comTemplate(specList,repSpecName)
%% Produce an ft.com template file in the templates directory

% This will be modified manually as well outside the program.
% Necessary paths also need to be reported to facilitate replacement:
% - data in
% - data out

%% Keep trying to get the file and run it; don't let the program proceed until cancel or success

            cancel = 0;
            success = 1;
            
            while ~(cancel==1 || success==0)
                
        %% Get/generate the file
    
            % Ask the user if they want to generate or pick an ft.com template

                response = menu('','Generate a basic ft.com file using basicFT1.com','Pick a template ft.com file','Cancel Processing');
                cd(specList.paths.templates)
                
                switch response
                    % From the proc_civm.com template file:
                    case 1
                            locs = [' ./representative_spectrum/',repSpecName,'/test.fid',... % data IN location
                                    ' ./representative_spectrum/',repSpecName,... % ft.com directory
                                    ' ./representative_spectrum/',repSpecName,... % data OUT location
                                    ' ',repSpecName];                      % spec name (used for naming files)
        %                     locs = [' ./test.fid',... % data IN location
        %                             ' ./',... % ft.com directory
        %                             ' ./']... % data OUT location
        %                    
                            success = system(['proc_civm.com',locs]);
                            file = [repSpecName,'_ft.com'];
                            path = pwd;
                            cancel = 0;

                     % Pick from existing template
                     case 2
                        %% Record location, make copy in proper location

                            % Get file info

                                [file,path] = uigetfile('*.com');
                                    dcfiles = dir([specList.paths.templates,'/*.com']);
                                    
                            % handle case where no file was selected
                                if file == 0
                                    success = 1;
                                    cancel = 0;
                                    continue % restart the loop
                                end
                                    
                            % Make sure the file (or copy) is in the templates
                            % directory

                                if ismember(file,{dcfiles.name}) % what if it goes by the same name??
                                    tempFile = 'externalFT.com';
                                    copyfile([path,file],[specList.paths.templates,'/',tempFile]) 
                                        %fileattrib([specList.paths.templates,'/',file],'+x','a');
                                else
                                    copyfile([path,file],[specList.paths.templates,'/',file])
                                    tempFile = file;
                                end
                                    fprintf(['\n\tCopying ''',file,''' to ',specList.paths.templates,' as ',tempFile,' to use as the template ft.com file.\n'])
                                    
                            % Tell user which file we're using

                               fprintf('\n\tUsing '),fprintf(tempFile),fprintf(' as the template ft.com file.\n')

                        % Modify the file to make it into the ft_com template file by removing spec-specific things

                            fdata = fileread([specList.paths.templates,'/',tempFile]);

                                replaceIn = './test.fid';
                                replaceOut = ['./',repSpecName','.ft'];

                                            fdata = regexprep(fdata,['(?<!#[^\n])',...  % not preceded by a # followed by any non-newline char (e.g. not commented out)
                                                                      '(?<=nmrPipe -in )',... % preceded by 'nmrPipe -in '
                                                                        '\S*',...       % expr is any number of consecutive non-whitespace chars
                                                                      '(?=\s\\)'],...   % followed by ' \'
                                                                     replaceIn);        % replacement for expr

                                            fdata = regexprep(fdata,['(?<!#[^\n])',... % not preceded by a # followed by any non-newline char (e.g. not commented out)
                                                                      '(?<=-out )',... % preceded by '-out '
                                                                        '\S*',...      % expr is any number of consecutive non-whitespace chars
                                                                      '(?= -ov)'],...  % followed by ' -ov'
                                                                     replaceOut);      % replacement for expr

                                % Write the file 

                                    path = specList.paths.rep_spec;
                                    newFile = [repSpecName,'_ft.com'];
                                    f = fopen([path,'/',newFile],'w'); % use the spec number as the filename
                                        fprintf(f,'%s',fdata);
                                    fclose(f);

                                % Overwrite file and path with new data (no
                                % need to reference to the old one)

                                    file = newFile;

                                % Make executable 

                                    fileattrib([path,'/',file],'+x','a');

                         cancel = 0;
                         
                    % Run the file 
                        cd(path)
                        success = system(file); % we want success = 0 (return value)

                 % Quit the program
                    case 3 
                        cancel = true;
                        success = 1;
                        return
                end                                               
                
            end
                                       

end