function fileInfo = getfidDotComFile(varargin)

% This function should accept an fid.com file or generate one from nmrPipe,
% then create individual fid files for each spectrum. Default is to run
% from the location of the file that will be processed.
%
% 'fromSpectrum',specNumber
% 'fromTemplate',pathToFile

%% Parse args

%     parms = {'fromSpectrum','fromTemplate'};
%     [ind] = strcmp(parms,varargin);
%
%     if ind(1)
%         fromSpectrum = varargin{:};
%     end

%% If generating a new file (default)

    % Using a specific spectrum (fromSpectrum):
            % go to the selected file
            
            % run bruk2pipe
                %!bruker
    
    % Using the first spectrum (default):
            % go to the first file on the list for that type
            
            % run bruk2pipe
                !csh -c 'bruker'
                
                % At this point, the user should Read Parameters and Save Script
        
    
%% fromTemplate

%                 templist = dir(pwd);
%
%
% % New section to check for procParm.txt
% %             % check for procParm.txt (new nmrPipe output)
% %                 if any(contains({templist.name},'procParm.txt'))
% %                     % Read Params
% %                 end
%
%
%             % If using the fid.com, extract the necessary parameters and store in a params sub struct
%                 if any(contains({templist.name},'fid.com'))
%                     pipepars = extractPipePars();
%                 else
%                     fprintf('\n\n\tNo fid.com file was found. Aborting. \n\n')
%                     return % kill the program
%                 end
      

    % Modify the fid.com file
        fidData = fileread('fid.com');
        fidData = regexprep(fidData,'sleep 5','sleep 0.001');
        
        % Write the new file to new scripts directory
            
            f = fopen('fid.com','w');
                fprintf(f,'%s',fidData);
            fclose(f);

        % Make it executable
            fileattrib('fid.com','+x','a');
            
        
    fileInfo.location = pwd;
    fileInfo.fullpath = [pwd,'/fid.com'];
    fileInfo.data = fidData;
end
