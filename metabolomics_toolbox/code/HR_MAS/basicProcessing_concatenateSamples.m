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
    ref = 1; % don't ref
    refppm = 0; % ppm
    refthresh = .004; % threshold for refspectra
    
    if ~isempty(varargin)
        ind = find(strcmp(varargin,'dataType'));
        if ~isempty(ind)
            dataType = varargin{ind+1};
        end
        
        ind = any(strcmp(varargin,'noRef')); % DON'T reference the spectra
        if ~isempty(ind)
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

    data = struct();
    
    for s = 1:length(studyInfo.sample)
        
        [~,ind] = ismember({dataType},{studyInfo.sample(s).expType.type});
            % NOTE: may want to allow passthrough of refSpec params
                data(s).spectra = HRMAS_nmr_runStdProc(studyInfo,s,ind,ref,refppm,refthresh);
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