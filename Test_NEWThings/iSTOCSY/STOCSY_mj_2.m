function [corr,covar,highCorrPPMs]=STOCSY_mj_2(target,X,ppm,threshold,intT,figureOption)

% STOCSY(target,X,ppm,threshold,figureOption,save)
% 
% Plots correlation/covariance projection of NMR spectrum to target
% chemical shift or response vector

% Arguments:
% 
% target       Chemical shift of driver peak or response vector with length equal to size(X,1)
% X            Data matrix of spectra
% ppm          Chemical shift vector corresponding to X
% threshold    Thresholding value (absolute value)
% figureOption Used to enable figures (e.g. for iSTOCSY; takes time; leave blank to turn off)
%                Generate figures: 'generateFigures'
%                Don't generate: 'false'
% save         Option for saving figures: 'save' (default, not used if
%               figureOption = 'false') 
%close all

%% Testing Params
%{
target = 8.435; 
X = XRALN; 
ppm = ppm; 
threshold = 0.85;
figureOption = 'generateFigures';

%}
save = 'save';
%% Parse and Check Input Params
tic
p = inputParser;
addRequired(p,'target',@(x) validateattributes(x,{'numeric'}, {'scalar'}));
addRequired(p,'X',@(x) validateattributes(x,{'numeric'}, {'2d'}));
addRequired(p,'ppm',@(x) validateattributes(x,{'numeric'}, {'row'}));

validPlotTypes = {'plot','stem'};
addParamValue(p, 'plotType', 'plot', @(x) any(validatestring(x,validPlotTypes)));

parse(p,target,X,ppm)

cellfun(@(f) evalin('caller',[f ' = p.Results.' f ';']), fieldnames(p.Results))

if length(target)==1
    [h,k]=min(abs(ppm-target));
    target_vect=X(:,k);
else
    target_vect=target;
end
    
if size(target_vect,1)~=size(X,1) && size(target_vect,2)==size(X,1)
    target_vect=target_vect';
end

%% Let us know what's running
fprintf(['Running STOCSY_mj on ' num2str(target) 'ppm...\n'])

%% Calculations 
corr=(zscore(target_vect')*zscore(X))./(size(X,1)-1);   
covar=(target_vect-mean(target_vect))'*(X-repmat(mean(X),size(X,1),1))./(size(X,1)-1);

%% Thresholding

thresh = threshold; % plus or minus
% {
%Option for Thresholding (block out if not wanted)
corrT = corr;
corrT(abs(corr)<thresh) = 0; %Set everything < |threshold| to zero -> corrT (thresholded correlation matrix) 
%}

%% Thresholding for intensity
%{
% Thresholding based on eyeballed baseline threshold intensity from mean
% spectrum (supplied as param):
%intT = 0.01;
corrT( (max(X,1) < intT)) = 0;

% Thresholding based on what was peakpicked?
%}
%% Figures
%%%%%% Block out the figure printing/saving to save tons of time in iSTOCSY (insert
%%%%%% space between % and {  to disable and turn on figures.
% {
if strcmp(figureOption, 'generateFigures')
        % set up colormap for raw STOCSY
        cmap=jet(100);
        lines_T=NaN(size(cmap,1),size(corr,2));
        lines_raw=lines_T;
        ind=1;
        for k=-1:2/size(cmap,1):.99
            lines_raw(ind,find(corr>k))=covar(find(corr>k)); % loop through the possible colormap values assigning data to them, sort of like thresholding for 100 different values with different colors. 
            ind=ind+1;
        end
        
        % set up colormap for thresh STOCSY
        cmap=jet(100);
        lines=NaN(size(cmap,1),size(corr,2));
        ind=1;
        for k=-1:2/size(cmap,1):.99
            lines_T(ind,find(corrT>k))=covar(find(corrT>k)); % loop through the possible colormap values assigning data to them, sort of like thresholding for 100 different values with different colors. 
            ind=ind+1;
        end
        %{
        % I never used the stem plot, it always made a blank figure for me
        if strcmp(plotType,'stem') 
        figure, stem(ppm,lines(1,:),'Color',cmap(1,:),'marker','none')
        hold on
        for k=2:size(cmap,1)
            stem(ppm,lines(k,:),'Color',cmap(k,:),'marker','none');
        end


        elseif
        %}       
        
        %% Generate the Figure (if opted)
        H = figure;
            if strcmp(plotType,'plot')
                ax(1) = subplot(3,2,[1 2]); % Plot the raw STOCSY data
                    %figure
                    plot(ppm,lines_raw(1,:),'Color',cmap(1,:))
                    hold on
                    for k=2:size(cmap,1)
                        plot(ppm,lines_raw(k,:),'Color',cmap(k,:));           
                    end
                    set(gca,'XDir','rev')
                    xlabel('Chemical Shift (ppm)')
                    if length(target)==1
                            ylabel(['Covariance with signal at ',num2str(target),' ppm'])
                        else
                            ylabel('Covariance with Y vector')
                    end
                    title(['STOCSY on Signal at ',num2str(target),' ppm'])
                ax(2) = subplot(3,2,[3 4]); % Plot the thresholded STOCSY data
                    %figure
                    plot(ppm,lines_T(1,:),'Color',cmap(1,:))
                    hold on
                    for k=2:size(cmap,1)
                        plot(ppm,lines_T(k,:),'Color',cmap(k,:));           
                    end
                    set(gca,'XDir','rev')
                    xlabel('Chemical Shift (ppm)')
                    if length(target)==1
                            ylabel(['Covariance with signal at ',num2str(target),' ppm'])
                        else
                            ylabel('Covariance with Y vector')
                    end
                    title(['Thresholded STOCSY (|R| > ' num2str(threshold) ') on Signal at ',num2str(target),' ppm'])     

                ax(3) = subplot(3,2,[5 6]); % Plot all the Spectra
                    %figure
                    plot(ppm,X,'b','LineWidth',0.01)  %%CHANGE HERE
                    set(gca,'XDir','reverse')
                    title('Aligned, Normalized Spectra for All Samples')
                    xlabel('Chemical Shift (ppm)')
                    ylabel('Relative Intensity')
                linkaxes(ax,'x')
            end
       %% Save the Figure (if opted and if figure is actually generated)
        thresholding = 'true'; % cannot produce the middle graph without thresholding
        if strcmp(save,'save') % only run if saving figures is desired
        filename = ['1DSTOCSY_on_',num2str(target),'ppm.fig'];
            if strcmp(thresholding,'true')
                filename = ['1DSTOCSY_on_',num2str(target),'ppm_threshold_',num2str(thresh),'.fig'];
            end
            savefig(H,filename,'compact')
            fprintf('Figure was saved...\n');
            close all
        else
            fprintf('Figures will not be saved...\n');
        end
end %

%% Get the interactions from the thresholded matrix
spectrum = mean(X,1); %% NOTE: the STOCSY graph will show the mean spectrum as the spectum, but this might well be changed to a covar of the dataset. 
% corrT is the thresholded matrix
%% Take the peaks from STOCSY output and get ppms in a separate file.
% the output file has the driver in column one, with the ppms of the
% responders in subsequenct columns

% Initialize Vars and Files
%interactions_outFile = ['STOCSY_cluster_for_',num2str(target),'ppm_threshold_',num2str(threshold),'.csv'];
%csvwrite(interactions_outFile,[]);
spec = zeros(size(spectrum));

%% Get the interactions from the thresholded correlation matrix

        [~,corrHighIndices,corrHighVals] = find(corrT); % this works. Take the list of indices of nonzero values and store in corrHighPeaks
    
        %Use the corrHighPeaks indices to find the corresponding values in spectrum matrix spec.
        spec(corrHighIndices) = spectrum(corrHighIndices); % put in spec only the peaks in the spectrum that the correlations would be projected onto.
        [~,peakMaxIndices] = findpeaks(abs(spec(:))); % take the max values of those peaks (responders) and store their indices, which line up with ppms
        highCorrPPMs = ppm(peakMaxIndices); % use the indices to get the ppms, print the responders next to the driver 
% NOTES: 
% - 'interactions refers to a previous iteration of this script. Here the
% passed result is called 'highCorrPPMs'.
% -The threshold must be high in order to get a reasonable number
% (<20) ppms responding. The algorithm basically:
% 1) takes the indices of the points on the spectrum (covar or mean) for which the correlation with the driver
% is above the specified threshold
% 2) lines them up with the corresponding values in the spectrum specified
% (covar or mean), like painting the regions of high enough correlation
% red, then takes the local maximum of each red region with regard to the
% corresponding spectrum intensity values. Then,
% 3) converts the indices of those local maxima to ppms and spits that list
% out in 'interactions' along with the driver (1st in the list)
% 
% I would like to develop this to choose the maxima of pre-picked peaks
% (i.e. from a peakpicking algorithm) if the correlated regions fall on one
% of those, as covariance may not be the best reflection of spectral intensity, which 
% databases like COLMAR rely on. The point is to be able to take the
% interactions output and directly enter into COLMAR or the like for
% putative compound ID. It may be possible eventually to iteratively search
% as a way of feeling out each peak's actual correlation threshold by
% responding dynamically to the number of results returned from COLMAR.
% Just a thought.

%%
        %}
        % %%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%
        % I don't see why this doesn't print a jet colorbar. Better just to
        % add manually on figure, if wanted. 
        %{
        t=colorbar;
        if length(target)==1
            set(get(t,'ylabel'),'String', ['Correlation with signal at ',num2str(target),' ppm']);
        else
            set(get(t,'ylabel'),'String', ['Correlation with Y vector']);
        end
        caxis([-1 1])
        %}
        toc
end
