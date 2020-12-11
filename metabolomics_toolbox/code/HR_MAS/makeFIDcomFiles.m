function [output] = makeFIDcomFiles(studyInfo,fidTemplate)
%%

    
     
%%        
    

    % (We'll need all of the files; no optimization takes place
    % until after the next step anyways)

    fnames = {studyInfo.sample(s).expType(t).files.name};

            fidcom2raw = fullPath2RelativePath(studyInfo.paths(s).fid_com,...
                       [studyInfo.paths(s).raw,'/',studyInfo.expTypes{s,t}],...
                       'useEscapeCharacters');

            fidcom2fid = fullPath2RelativePath(studyInfo.paths(s).fid_com,...
                       studyInfo.paths(s).fid,...
                       'useEscapeCharacters');

    cd(studyInfo.paths(s).fid_com)

    for f = 1:length(fnames)
        fidcomname = make_fidDotCom_customPaths(fileInfo.data,...
                                                            fidcom2raw,...
                                                            fidcom2fid,...
                                                            studyInfo.paths(s).fid_com,...
                                                            fnames{f});
        system(fidcomname);
    end
            
%% Outputs

    output.commands = fidcomname;
    output.files = dir();
            
end