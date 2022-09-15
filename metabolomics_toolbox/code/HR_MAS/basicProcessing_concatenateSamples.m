function [catData] = basicProcessing_concatenateSamples(studyInfo,varargin)
%%
% 
%     Do basic processing on 1D CIVM data, including concatenation into a structure.
%     studyInfo 
%       Varargin:
%                       
% 
%% Handle varargin
    %     NOTE: may want to allow passthrough of refSpec params
    %       

    %     sampOrder = {'civm_ncrassa_qax_02_media',...
%                 'civm_ncrassa_qax_02_media_qa',...
%                 'civm_ncrassa_qax_02_media_org',...
%                 'civm_ncrassa_qax_02_media_org_qa_run1',...
%                 'civm_ncrassa_qax_02_media_org_post_qa_spike_iconRun'
%                 };


%% Do basic processing on everything
    
    dataType = 'noesypr1d'; % typical CIVM 1d expt
    ref = 1; % ref
    refppm = 0; % ppm
    refthresh = .004; % threshold for refspectra
    maxWithin = []; % null initiation; maxWithin to pass ref window to ref function
    
    if ~isempty(varargin)
        ind = find(strcmp(varargin,'dataType'));
        if ~isempty(ind)
            dataType = varargin{ind+1};
        end
        
        ind = any(strcmp(varargin,'noRef')); % DON'T reference the spectra
        if ind
            ref = 0;
        end
        
        ind = find(strcmp(varargin,'refppm'));
        if ~isempty(ind)
            refppm = varargin{ind+1};
        end
        
        ind = find(strcmp(varargin,'refthresh'));
        if ~isempty(ind)
            refthresh = varargin{ind+1};
        end
        
        ind = find(strcmp(varargin,'maxWithin'));
        if ~isempty(ind)
            maxWithin = varargin{ind+1};
        end
        
        
%         ind = find(strcmp(varargin,'refThreshold'));
%         if ~isempty(ind)
%             noRef = 0;
%         end
%         
%         ind = find(strcmp(varargin,'refOffset'));
%         if ~isempty(ind)
%             noRef = 0;
%         end
%         
%         ind = find(strcmp(varargin,'refManualPick'));
%         if ~isempty(ind)
%             noRef = 0;
%         end
    end

    %%
    data = struct();
    
    for s = 1:length(studyInfo.sample)
        
        [~,ind] = ismember({dataType},{studyInfo.sample(s).expType.type});
        if ind == 0
            warning(['Supplied or default argument for ''dataType'' "',dataType,...
                    '" was not found in current sample "', studyInfo.sample(s).name,'"',...
                    ' ... Trying data type found ("',studyInfo.sample(s).expType(1).type,'")....',...
                    ' If this works, consider changing the call to basicProcessing_concatenateSamples()',...
                    ' to include name-value pair argument: ''dataType'', ''',studyInfo.sample(s).expType(1).type,''''])
                dataType = studyInfo.sample(s).expType(1).type;
                [~,ind] = ismember({dataType},{studyInfo.sample(s).expType.type});
                if ind == 0
                    error(['dataType ',studyInfo.sample(s).expType(1).type,' not found. Quitting.'])
                end
        end
            % NOTE: may want to allow passthrough of refSpec params
            
%                 data(s).spectra = HRMAS_nmr_runStdProc(studyInfo,s,ind,ref,refppm,refthresh,maxWithin);
                data(s).spectra = HRMAS_nmr_runStdProc(studyInfo,s,ind,...
                                                        'doRef',ref,...
                                                        'refppm',refppm,...
                                                        'refthresh',refthresh,...
                                                        'maxWithin',maxWithin);
                                                    
    end
    
    % Unlist the data struct one step
        
        sspectra = vertcat(data(:).spectra);
                
%% Concatenate the data

    % Get all the spectra in the same struct level
    
        % Unlist the data struct again (not sure how else to do this)
        catData = sspectra(1).spectra;
        
        if length(studyInfo.sample)>1
            for s = 2:length(studyInfo.sample)
                catData = [catData,sspectra(s).spectra];
            end
        else
            catData = sspectra(s).spectra;
        end

%         catData(~ismember({catData.experiment},{dataType})) = []; % remove empties
        catData(cellfun(@isempty,{catData.real})) = []; % remove empties
        
        [~,inds] = sort([catData.startTime]); % sort by time
            catData = catData(inds);    
        cd(studyInfo.sample(1).paths.sample), cd ..
    
end