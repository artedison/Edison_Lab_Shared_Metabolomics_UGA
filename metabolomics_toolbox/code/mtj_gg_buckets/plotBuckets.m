function plotBuckets(matrix,ppm,bins,varargin)

%% plotBuckets

    % Author: MTJ
    % Version: 0.2
    % Tested on Matlab Version R2020a
    % Date: JUL2020
    %
    % Description:
    %       
    %       Generates a overlayed spectral plot with buckets drawn. Options 
    %       included for expandedBuckets (only outlines) and plotting in a 
    %       subplot context. 
    %
    % Input:
    %       matrix:     spectral matrix
    %       ppm:        ppm vector
    %       bins:       array of bucket boundaries (ppm values, 2 x n)
    %       varargin:     
    %                   'expandedBuckets' - use for plotting expanded
    %                                       buckets (see
    %                                       optOB_expandBuckets() and 
    %                                       expandBucketBounds() ). This
    %                                       causes buckets to be
    %                                       transparent color with black
    %                                       borders. 
    %                   'subplot'         - if plotting into a subplot (as
    %                                       in plotOptBucket_optResult() )
    %                   'title'           - optional title string
    %                                       (name-value pair)
    %                   
    % Output:
    %       
    %       Plot with spectra and buckets drawn on them. 
    %
    % Log:
    %
    % Example run:
    %           
    %       plotBuckets(X,ppm,[1.2,3.4,7.2;1.24,4,7.3],'expandedBuckets','title','Plot Title')
    %


%% Handle varargin
    
    expandedBuckets = false;  % default
    subplot = false;
    titlestr = 'Buckets';
    
    if ~isempty(varargin)
        if ismember('expandedBuckets',varargin)                          
            expandedBuckets = true;
        end
        if ismember('subplot',varargin)                          
            subplot = true;
        end
        if ismember('title',varargin)      
            [~,ind] = ismember('title',varargin);      
            titlestr = varargin{ind+1};
        end     
    end    
    
%% Do the Plotting
    
    % If you want to make a new figure (vs. putting in an existing axes)
    if ~subplot
        figure, 
    end
        hold on
            
            plot(ppm,matrix)    
                %scatter(optOB_out.results(n).matchingPeaks,max(optOB_out.results(n).matchingPeakInts ,[],1),'o','k')
            set(gca,'XDir','reverse')
            xlabel('Chemical Shift (ppm)')
            ylabel('Signal Intensity')
            %set(gca,'xlim',xlims)
            %set(gca,'ylim',ylims)
            %set(gca, 'YTickLabel',[])
            title(titlestr)

    %% Highlight the bins (changes if expanded)

        if expandedBuckets
            highlightROIs( bins', max(matrix(:)) ,'color','none','edgeColor','k')
        else
            highlightROIs( bins', max(matrix(:)) ,'color','r')
        end
        
        hold off    
        
end