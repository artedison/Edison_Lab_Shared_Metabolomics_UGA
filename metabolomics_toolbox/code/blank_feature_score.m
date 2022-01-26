function      [per_spectra_ppick_output] = blank_feature_score (matrix, blank_matrix, ppm, ppick_output,  per_spectra_ppick_output, sorted_bucket_list)


% This function is part of the package of functions developed by MJ and GG. 
% should be used only after bins and peak peaked functions
% It is only for deteermining peaks that are in bblanks and compare against the blank remooved dataset 


%% see at the end for sequeential fucntions and usage%




%% Determines how many Peaks from blanks are present in each bin 
% 0 - no peaks detected in each bucket >0 - one or more peaks in each
% bucket

for ii=1:size (sorted_bucket_list,1)
    
    low_buck = sorted_bucket_list(ii,1);
    high_buck =sorted_bucket_list(ii,2);
    
    for i = 1:size(ppick_output.blank_shifts,2)
        
        per_spectra_ppick_output.blanks.blanks_flag(ii,i) =...
            any(...
            ppick_output.blank_shifts(:,i) >=low_buck...
            &...
            ppick_output.blank_shifts(:,i)<=high_buck);
        
        
    end
end

per_spectra_ppick_output.blanks.blanks_flag = sum( per_spectra_ppick_output.blanks.blanks_flag,2);


%% Getting the max blank peak for each bucket and its respective ppm value. If no peak detected the mean ppm for each peak is assigned

for i = 1:size(blank_matrix,1)
    for    ii=1:length (sorted_bucket_list)
        
        if  isempty(...
                blank_matrix(i,...
                ppick_output.blank_shifts>=sorted_bucket_list(ii,1)...
                &...
                ppick_output.blank_shifts<=sorted_bucket_list(ii,2) ) )  == 1
            
            per_spectra_ppick_output.blanks.blank_ints(i,ii)= nan;
            
            per_spectra_ppick_output.blanks.blank_shift(1,ii)= per_spectra_ppick_output.max_iteration.gap_f_median_ppm(ii);
            
        else
            [ per_spectra_ppick_output.blanks.blank_ints(i,ii),pp] = ...
                max(...
                    ppick_output.blank_ints(i,...
                    ppick_output.blank_shifts>=sorted_bucket_list(ii,1)...
                    &...
                    ppick_output.blank_shifts<=sorted_bucket_list(ii,2) )...
                ,[],2) ;
            
            temp_ppm =  ppick_output.blank_shifts(:,...
                ppick_output.blank_shifts>=sorted_bucket_list(ii,1)...
                &...
                ppick_output.blank_shifts<=sorted_bucket_list(ii,2) );
            
            per_spectra_ppick_output.blanks.blank_shift(1,ii) = temp_ppm(:,pp);
            
        end
    end
end


%%


per_spectra_ppick_output.blanks.blank_ints_flag =...
    mean( per_spectra_ppick_output.blanks.blank_ints)...
    >=...
    mean(per_spectra_ppick_output.max_iteration.gap_f_ints)./3;

%%
colors = [69 69 69 ;... %samples
                111 134 149;...%bucket boundaries
                0, 0, 0;... %1/3 mean
                242 112 149;... %blanks
                227 86 144;... %max blank
                54 129  127;... % peak picked blank max
                ]./255;

figure, hold




plotr(...
    ppm,...
    matrix,...
    '--', 'color', colors(1,:))

line( [sorted_bucket_list(:) sorted_bucket_list(:)],  get(gca, 'ylim'),'Color',colors(2,:), 'LineWidth',0.75 ) % boundaries

plotr(...
    ppm,...
    mean(matrix)./3,...
    'color', colors(3,:), 'LineWidth', 3)

plotr(...
    ppm,...
    blank_matrix,...
    '-', 'color', colors(4,:))

plotr(...
   ppm,...
    mean(blank_matrix),...
    'color', colors(5,:), 'LineWidth', 3)

plotr (...
    per_spectra_ppick_output.blanks.blank_shift,...
    mean(per_spectra_ppick_output.blanks.blank_ints),...
    '.','MarkerSize',16,  'color', colors(6,:))


text(...
    per_spectra_ppick_output.max_iteration.gap_f_median_ppm,...
    mean(per_spectra_ppick_output.blanks.blank_ints).*0-mean(matrix(:)/0.05),...
    string(double(per_spectra_ppick_output.blanks.blanks_flag)),...
    'color','red','clipping','on')

text(...
    per_spectra_ppick_output.max_iteration.gap_f_median_ppm,...
    mean(per_spectra_ppick_output.blanks.blank_ints).*0-mean(blank_matrix(:)/.025),...
    string(double(per_spectra_ppick_output.blanks.blank_ints_flag)),...
    'color','blue','clipping','on')


addReasonableLegend(...
    [{'Samples'},{'Bucket Boundaries'}, {' Blank threshold (sample mean / 2) '}, {'Blanks'}, {'Max of all blanks'}, {'Peak picked blank features'}],...
    colors, 'addBox','yes', 'textSize', 9,'lineWidths', 5)

   %% Optimize Peak Picking (threshold; for representative spectrum)
%         
%                                 matrix = wrk_data.XRA;
%                                 ppm = wrk_data.ppmR;
%                                      optimize_Peakpick1D(matrix,ppm,'var',0.05:0.05:0.25,'Complex');   
%                             %% Do the actual Peak Picking (for representative spectrum)
%                                      peaks = struct();
%                                      [peaks.ints, peaks.shifts]= Peakpick1D(matrix ,ppm,'max',0.2,'Complex');
%                        
%                             %% Automatic Binning (Bucketing)
%  
%                             %% Generate buckets using a range of both params
% 
%                                 sb = 0.002:0.002:0.008;
%                                 sl = 0.3:0.1:0.6;
%                                 [optOB_out] = optimize_optBucket(matrix,ppm,sb,sl);
%                             %% Filter out the bins with no peaks
%                                [optOB_out] = filterBuckets_Peaks_opt(ppm,optOB_out, peaks);        
%                             %% Plot the results of optOB
%                                  [optOB_out] =  plotOptBucket_optResult(matrix,ppm,optOB_out,[3.6584    4.0], [min(matrix(:)) max( matrix(:,(ppm>3.6584  & ppm<4.0)), [],'all' ) ]);
%                          %% Peak pick every spectra independently, so that each spectrum contributes to the peak shape between the boundaries
%                                 peakthresh = 0.15 % adjust for moore or fewer peaks
%                                 mode = 'Complex'      
%                                 [itpeaks] = Peakpick1D_per_spectra(matrix,ppm,peakthresh,mode)
%                          %%  Expand buckets to each of the bins boundaries
%                                  [optOB_out] = expandBucketBounds (optOB_out, matrix, ppm, 'plotResult');
%                         %%  Manual Refinement of the boundaries
%                             %% [optOB_out] = refineBuckets(matrix,ppm,optOB_out,5);    
%                                  [optOB_out] = refineBuckets_GG2(matrix,ppm,optOB_out,itpeaks,peaks, 'expandedBuckets');    
%                              %% saving variables of the workspace 
%                                save('post_buckets_refined_xx_xx_2022.mat')
% %% This calculates the number of peaks in each bin in the original Peakpick1D_per_spectra data
%         % Calculates the max peak within each bin and its chemical shift 
%         % gap fills both ppm and intensities of peaks that were no present
%         % in the bin or not picked (not detected/below baseline)
%  matrix = wrk_data.XRA;
%  ppm = wrk_data.ppmR; 
%  buckets = optOB_out.refinedBuckets.refinedBuckets;
%  ppeak_struct = itpeaks;
% [itpeaks] = peaks_per_bin (matrix, ppm,buckets, ppeak_struct);
% %% Align metric
%  buckets = itpeaks.max_iteration.sorted_bucket_list;
%  ppeak_struct = itpeaks;
%  metadata = Td_pd_ugt;
% itpeaks = align_metric ( matrix, ppm, ppeak_struct, buckets, metadata );
% clearvars matrix ppm 
% %% getting blank matrix peak picked
% blank_matrix = wrk_data.X_blks;
% ppm = wrk_data.ppm; 
% [peaks.blank_ints, peaks.blank_shifts]= Peakpick1D(blank_matrix ,ppm,'max',0.9,'Complex');
% %% scoring blank peaks 
% spectra_struct = wrk_data.X_pd_ugt;
% pick_output = peaks;
% per_spectra_pick_output = itpeaks;
% sorted_bucket_list = itpeaks.max_iteration.sorted_bucket_list;
%  [itpeaks] = blank_feature_score (spectra_struct, blank_matrix, ppm, pick_output,  ...
%                                                         per_spectra_pick_output, sorted_bucket_list)
% clearvars spectra_struct ppick_output per_spectra_pick_output sorted_bucket_list
% clearvars blank_matrix ppm

