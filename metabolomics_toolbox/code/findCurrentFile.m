function [parentpath,filepath,filename] = findCurrentFile()
% findCurrentFile()
%
% Description: Finds and returns the local parent directory (one up from
% file location), file location, and filename for the workflow being run
% within the MATLAB editor. Intended use is for building directory tree
% structures of local filepaths at the beginning of a workflow run to avoid
% dependence on manual input for adding files to path and navigation.  
% 
% modified from: https://www.mathworks.com/matlabcentral/answers/81148-get-path-from-running-script#answer_395208
% related: matlab.desktop.editor.getAll from https://www.mathworks.com/matlabcentral/answers/119438-how-can-i-create-a-list-of-m-files-currently-open-in-the-editor#answer_126434
%
% Useage (called from active workflow):
%         [paths.project,paths.scripts,paths.thisWorkflow] = findCurrentFile(); 
%             paths.results = [paths.project,'/','results'];
%             paths.data = [paths.project,'/','data'];
%                 paths.raw = [paths.data,'/','raw'];
%                 paths.processed = [paths.data,'/','processed'];
%             paths.metadata = [paths.project,'/','metadata'];
%
% MTJ 2020

        wd = matlab.desktop.editor.getActiveFilename;
        spltwd = strsplit(wd,'/'); % navigate to the folder the workflow sits in
        
        filepath = join(spltwd(1:end-1),'/');
            filepath = filepath{:};
        parentpath = join(spltwd(1:end-2),'/');
            parentpath = parentpath{:};
        filename = spltwd{end};
        
end
