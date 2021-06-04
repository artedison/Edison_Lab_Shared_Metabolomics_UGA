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
    
%     dataType = 'noesypr1d';
% 
%     if ~isempty(varargin)
%         ind = strcmp(varargin,'dataType');
%         if ~isempty(ind)
%             dataType = varargin{ind+1};
%         end
%     end

    data = struct();
    
    for s = 1:length(studyInfo.sample)
        
        [~,ind] = ismember({'noesypr1d'},{studyInfo.sample(s).expType.type});
            % NOTE: may want to allow passthrough of refSpec params
                data(s).spectra = HRMAS_nmr_runStdProc(studyInfo,s,ind);
                
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

        [~,inds] = sort([catData.startTime]);
        catData = catData(inds);
        cd(studyInfo.sample(1).paths.sample), cd ..
    
end