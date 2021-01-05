function [output] = makeFIDcomFiles(specList,fidTemplate)
%%

    
     
%%        
    
    % (We'll need all of the files; no optimization takes place
    % until after the next step anyways)

    fnames = {specList.files.name};

            fidcom2raw = fullPath2RelativePath(specList.paths.fid_com,...
                       specList.paths.raw,...
                       'useEscapeCharacters');

            fidcom2fid = fullPath2RelativePath(specList.paths.fid_com,...
                       specList.paths.fid,...
                       'useEscapeCharacters');

    cd(specList.paths.fid_com)
    
    for f = 1:length(fnames)
        fidcomname = make_fidDotCom_customPaths(fidTemplate.data,...
                                                            fidcom2raw,...
                                                            fidcom2fid,...
                                                            specList.paths.fid_com,...
                                                            fnames{f});
    end
            
    %% Run the .com files
            temp2fidcom = fullPath2RelativePath(specList.paths.templates,...
                       specList.paths.fid_com,...
                       'useEscapeCharacters');
        
        cd(specList.paths.templates)
        system(['runFIDfiles.com ',temp2fidcom]);
    
%% Outputs

    output.commands = fidcomname;
    output.files = dir();
            
end