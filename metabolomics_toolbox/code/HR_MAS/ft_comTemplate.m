function [file,path,cancel] = ft_comTemplate(specList,repSpecName)
%% Produce an ft.com template file in the templates directory

% This will be modified manually as well outside the program.
% Necessary paths also need to be reported to facilitate replacement:
% - data in
% - data out

%% Ask the user if they want to generate or pick an ft.com template

            response = menu('','Generate a basic ft.com file using basicFT1.com','Pick a template ft.com file','Cancel Processing');
                
                            cd(specList.paths.templates)

%% Get/generate the file

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
                        system(['proc_civm.com',locs])
                        file = [repSpecName,'_ft.com'];
                        path = pwd;
                        cancel = 0;
                        
                 % Pick from existing template
                 case 2
                    % Record location, make copy in proper location
                        % Alert user of this fact

                        [file,path] = uigetfile('*.com');
                            %dcfiles = dir([specList.paths.templates,'/*.com']);
                        
                        % Make sure the file (or copy) is in the templates
                        % directory
                        if ~ismember(file,{dcfiles.name})
                            copyfile([path,file],specList.paths.templates)
                            fprintf(['\n\tCopying ''',file,''' to ',specList.paths.templates,' to use as the template ft.com file.\n'])
                        else
                            fprintf(['\n\tUsing ''',file,''' as the template ft.com file.\n'])
                        end
                        
                        % Modify the file to make it into the ft_com
                        % template file by removing spec-specific things
                        
                        fdata = fileread([specList.paths.templates,'/',file]);
                        fdata = regexprep();
                        fdata = regexprep();
                        fdata = regexprep();
                        
                        cd(path)
                        system(file)    
                case 3 
                    cancel = true;
                    return
            end 

end