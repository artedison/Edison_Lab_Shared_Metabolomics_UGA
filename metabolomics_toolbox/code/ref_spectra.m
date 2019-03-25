function [spectra, offsetppm] = ref_spectra(spectra,thresh,offsetppm,varargin)

    % Author: Edison Lab
    % Version: 0.2
    % Tested on Matlab Version R2017b
    % Date: 25MAR2019
    %
    % Description:
    %       This function provides options for referencing NMR spectra
    %       (e.g. to DSS or another chemical shift standard peak). NOTE:
    %       referencing spectra requires interpolation of new ppm vectors.
    %       This can result in slight variations in maximum peak height for
    %       a peak, typically 1 ppm resolution unit difference in either
    %       direction. 
    %
    % Input:
    %   Required
    %       spectra : structure output from Setup1D containing spectral
    %                   data and ppm vectors for each spectrum
    %   Optional
    %       thresh : peak picking threshold. Default: 0.01.
    %       offsetppm : Expected location of reference peak. Pick the peak 
    %                   in each spectra closest to offsetppm. Providing 
    %                   this suppresses manual picking of reference peak. 
    %       varargin 
    %       	'testThreshold' : providing this argument allows the result
    %                               of peakpicking to be assessed in order
    %                               to optimize thresh without having to
    %                               pick any peaks/do the referencing
    %                               calculation.
    %
    % Output:
    %       spectra : the referenced spectra and adjusted ppm vectors
    %       offsetppm: this facilitates iterative running for optimization
    %                   purposes, as well as allows for recordkeeping.
    %
    % Log:
    %       Ver 0.0: old toolbox function
    %       Ver 0.1: MJ and YW edit 10/10/2018
    %       Ver 0.2: testThreshold option added MJ FEB2019
    %
    % Example run:
    %
    
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
                figure, hold on
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