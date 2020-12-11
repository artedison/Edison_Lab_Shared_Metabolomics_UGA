function [output] = makeAndRunFT_comFiles(genFilename,runFilename,studyInfo)
%%


%% Go to the template directory (we'll run the scripts from here)

        cd(studyInfo.paths_sample(s).pipe_templates)
                    
    %% Generate relative filepaths   
 
        temps2fid = fullPath2RelativePath(studyInfo.paths_sample(s).pipe_templates,...
            studyInfo.paths(s).fid,...
            'useEscapeCharacters');                    

        fid2temp = fullPath2RelativePath(studyInfo.paths(s).fid,...
            [studyInfo.paths_sample(s).pipe_templates,'/template_ft.com'],...
            'useEscapeCharacters'); 

        fid2ftcom = fullPath2RelativePath(studyInfo.paths(s).fid,...
            studyInfo.paths(s).ft_com,...
            'useEscapeCharacters');

        temps2ftcom = fullPath2RelativePath(studyInfo.paths_sample(s).pipe_templates,...
            studyInfo.paths(s).ft_com,...
            'useEscapeCharacters');
                    
        temps2ft = fullPath2RelativePath(studyInfo.paths_sample(s).pipe_templates,...
            studyInfo.paths(s).ft,...
            'useEscapeCharacters');
        
    %% Generate the individual ft.com files

                       cmd = [genFilename,...
                                ' ',temps2fid,...   % directory containing the fid files
                                ' ','specNumber',...  % str to replace in input file
                                ' ',fid2temp,...   % input file relative path
                                ' ',fid2ftcom];      % output files relative path (ft_com directory)

                       system(cmd)

    %% Run the ft.com files

                       system([runFilename,' ',temps2ftcom]);
                       
    %% Output
        
            output.commands.generate = cmd;
            output.commands.run = [runFilename,' ',temps2ftcom];
            output.relativePaths.temps2fid = temps2fid;
            output.relativePaths.temps2fid = temps2ft;
            output.relativePaths.temps2fid = fid2temp;
            output.relativePaths.temps2fid = fid2ftcom;
            output.relativePaths.temps2fid = temps2ftcom;
            output.files.ft_com = dir(temps2ftcom);
            output.files.ft = dir(temps2ft);
            
end