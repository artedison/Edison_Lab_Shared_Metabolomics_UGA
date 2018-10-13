function [commonFiles,filesToCopy] = mergeDirectories(dependencies,localpath,destination)

%% Set up params, etc.
%     localpath = '/Users/mjudge/';
%     destination = '/Users/mjudge/Edison_Lab_Shared_Metabolomics_UGA/metabolomics_toolbox/code/HR_MAS';
%     dependencies = resfilelist;

% %     flist=matlab.codetools.requiredFilesAndProducts({'STEP_1_processing_combine_samples.m','STEP_2_ridge_tracing.m','STEP_3_combining_ridges.m','STEP_4_plotting.m'});
% %     dependencies={};
% %     for i=1:length(flist)
% %       if ~any(regexp(flist{i},'.mat'))
% %         dependencies{end+1}=strrep(flist{i},'/home/yuewu/Documents/github/github.code/clone','');
% %       end
% %     end
% %     dependencies{end+1}='/Edison_lab_UGA/metabolomics_toolbox/Test_NEWThings/yuewu/filesys_cons.m';    
    
%% Compare dependencies to the existing files in the specified repository
    existingFiles = dir('**/*.m');

    % % Concatenate the path with the filename for each file
    %     s1 = {existingFiles.folder};
    %     s2 = {existingFiles.name};
    %     s3 = repmat({'/'},1,length(s1));
    % 
    %     files = cellfun(@(s1,s2,s3)[s1,s3,s2],s1,s2,s3,'UniformOutput',false);
    % Parse dependency names
        existingNames = {existingFiles.name};
        splitPath = cellfun(@(s) strsplit(s, '/'), dependencies, 'UniformOutput', false);  
        dependencyNames = cellfun(@(s) s{end}, splitPath, 'UniformOutput', false);
        
    % Intersect the lists to see which files exist so we can ignore them
        [commonFiles,~,ignoreInds] = intersect(existingNames,dependencyNames);
        
    % Report
        fprintf(['\n\n',...
                'The following functions already exist in the provided repository and will be ignored:',...
                '\n\n'])
        disp(commonFiles);

%% Copy the files over

    for i = setdiff(1:length(dependencyNames),ignoreInds)
        copyfile( [localpath,dependencies{i}], destination );
    end

    % Report 
        filesToCopy = dependencyNames(setdiff(1:length(dependencyNames),ignoreInds));
        fprintf(['\n\n',...
                'The following functions were added to the specified repository:',...
                '\n\n'])
        disp(filesToCopy);