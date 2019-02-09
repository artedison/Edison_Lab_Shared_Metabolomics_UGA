function [spectra, offsetppm] = ref_spectra(spectra,thresh,offsetppm,varargin)
% ************
% Reference spectra.
% Required arguments
%     spectra - Structure array containing fields: real, ppm
% Optional arguments:
%     thresh - 0.01 (default)
%     offsetppm - If set, don't manually pick the reference from the first spectra.
%     Instead, pick the peak in each spectra closest to offsetppm.
% OUTPUT: spectra and offsetppm of referenced spectral
%% it is an initial function from the toolbox
%% MJ and YW edit it 10/10/2018
% *******************
%% Check to see if we are testing the threshold
    if any(strcmp(varargin,'testThreshold'))
        testThreshold = true;
    else
        % do nothing
    end
    
%% Deal arguments
if ~exist('thresh','var')
    thresh = 0.01;
end

%% Test to see what your threshold should be
    if exist('testThreshold','var')
        % MJ 3JAN2018: Run this section to see if your threshold is low enough
        % (there should be a red circle on all of your Ref peaks; if not, lower
        % the threshold):
        peaks = cell(size(spectra,2),1);
        for i=1:length(spectra)
            peaks{i}=peakpick(spectra(i).real,30,thresh*max(spectra(i).real));
        end

        figure;
        hold on
        cmap=colormap(lines(length(spectra)));


        for i=1:length(spectra)
            plot(spectra(i).ppm,spectra(i).real,'Color',cmap(i,:))
            plot(spectra(i).ppm(peaks{i}(:)),spectra(i).real(peaks{i}(:)),'ro')
        end
        currentppm = spectra(1).ppm;
        % Set the limits to +/- 0.5 ppm from the offset (ref peak), but make
            % sure those points actually exist in the ppm vector).
            set(gca,'xlim',currentppm(matchPPMs([offsetppm - 0.5,offsetppm + 0.5],currentppm)));
            title('Peaks picked during ref_spectra (test mode)','Interpreter','none')
            
        msgbox(['This algorithm requires the reference peak to be peak-picked in all spectra. ',...
                'If too many peaks are picked (red circles), then raise the threshold. ',...
                'If the reference peak is not picked in all spectra, lower the threshold. ',...
                'Ideally, noise is not selected, but this is allowed. ',...
                'Current Threshold: ',num2str(thresh),' . ',...
                'Reference Peak (provided): ',num2str(offsetppm),' ppm.'])
    else
%% If we are not testing, complete the referencing:
        %% Manually pick reference peak and do the referencing
            if ~exist('offsetppm','var') || isempty('offsetppm')
                i=1;
                peaks=peakpick(spectra(i).real,30,thresh*max(spectra(i).real));
                plot(spectra(i).ppm,spectra(i).real,'k'), hold,
                plot(spectra(i).ppm(peaks),spectra(i).real(peaks),'ro')
                set(gca,'xlim',[-.7,.7])
                set(gca,'xdir','rev')
                [x,y]=ginput(1);
                close
                [~,idx]=min(abs(x-spectra(i).ppm(peaks)));
                TSP_ppm(i)=spectra(i).ppm(peaks(idx));
                spectra(i).ppm=spectra(i).ppm-TSP_ppm(i);
                offsetppm = TSP_ppm(i);
            else
                peaks=peakpick(spectra(1).real,30,thresh*max(spectra(1).real));
                [~,idx]=min(abs(offsetppm-spectra(1).ppm(peaks)));
                TSP_ppm(1)=spectra(1).ppm(peaks(idx));
                spectra(1).ppm=spectra(1).ppm-TSP_ppm(1);
            end
        %% Pick all peaks

            for i=2:length(spectra)
                spectra(i).ppm=spectra(i).ppm-TSP_ppm(1);
                peaks=peakpick(spectra(i).real,30,thresh*max(spectra(i).real));
                [~, idx] = min(abs(spectra(i).ppm(peaks)));
                TSP_ppm(i)=spectra(i).ppm(peaks(idx));
                spectra(i).ppm=spectra(i).ppm-TSP_ppm(i);
            end

            figure;
            hold on
            cmap=colormap(lines(length(spectra)));


            for i=1:length(spectra)
                plot(spectra(i).ppm,spectra(i).real,'Color',cmap(i,:))
            end
            set(gca,'XDir','reverse')
    end

end