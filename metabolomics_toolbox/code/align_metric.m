function [ppeak_struct ] = align_metric ( matrix, ppm, ppeak_struct, bucket_list, metadata);

%% alignment scale 
% This function generates 3 plots and a vector of the euclidean distances
% of each peak in one given bucket, across the whole spectra.
%
% Function requirement:
%       Peakpick1D_per_spectra - This generates the structure ppeak_struct
%       peaks_per_bin - this output into the structure ppeak_struct
%       refineBuckets_GG - this creates the refined bucket list 


% Inputs:
%    matrix  - -  where rows are samples and columns are ppms
%   ppm      - -  ppm vector with the same number of columns as matrix and
%                                    one row
%   buckets - - A matrix of 2 columns where each row specifies one bucket
%   ppeak_struct.max_iteration - containing - ppeak_struct.max_iteration.sorted_bucket_list, 
%                                                                  ppeak_struct.max_iteration.ints, 
%                                                                  ppeak_struct.max_iteration.shift,
%                                                                  ppeak_struct.max_iteration.max_elements_per_bin, 
%                                                                  ppeak_struct.max_iteration.gap_f_shift, 
%                                                                  ppeak_struct.max_iteration.gap_f_ints
%
% Outputs:
% Fig1 - 'Alignment scale per bin and picked peaks' - This plots all the
% spectra in black, each peak picked, each max peak within each bucket and
% gap filled peaks. The text in the plot represent the alignment value for
% each bin
%
% Fig2 - 'Alignment metric summary' - a representation of the alignment scale across all the bins 
%
%
% Fig3 - 'Defining a quality of alignment to scale values ranges' - creates
% ranges of the alignment scale and plots the prevalence of each of the
% ranges across all of the buckets while associating a quality to each
% range
%
% Note: Structures names can be user defined, however substructures and
% fields need to remain as defined by each function and dependency
% functions


%% Calculating euclidean distances from each point to the median of the data within a bin 

for i = 1:size( bucket_list,1)
    
    matx =ppeak_struct.max_iteration.shifts(ppeak_struct.max_iteration.shifts(:,i)~=0,i); % only of the max of the peak picked - does not take into acocunt the gap filled peaks
     
    
    
    euc_dist = ...
        pdist2(...
        matx,...
        mean(matx),...
        'euclidean' );
    
    
    ppeak_struct.align_scale.mean_dif_mean (:,i)=mean(euc_dist, 1); % mean distance per bucket
    ppeak_struct.align_scale.var_dif_mean (:,i) = std(euc_dist,1); % standard deviation per bucket
  
    clear euc_dist maxN
end




%% adds the std to each mean so it increases the differences between bins (more variance = higher values)

ppeak_struct.align_scale.metric =round( (ppeak_struct.align_scale.mean_dif_mean...
                                                                    ).*10000,...
                                                                    2,'decimals');


%%

figure,
hold
plotr(ppm, matrix,'k')

for i= 1:size(matrix,1)
    plotr( ppeak_struct.max_iteration.gap_f_shifts(i,:),ppeak_struct.max_iteration.gap_f_ints(i,:),'.', 'Color', [51, 108, 215]./255,'MarkerSize',30)
    hold on
end

line( [bucket_list(:) bucket_list(:)],  get(gca, 'ylim'),'Color', [111, 134, 149]./255, 'LineWidth',1 ) % boundaries
hold on

for i= 1:size(matrix,1)
    plotr( ppeak_struct.iteration(i).shifts,  ppeak_struct.iteration(i).ints,'o', 'Color', [242, 112, 89]./255,'MarkerSize',9)
    hold on
end

for i= 1:size(matrix,1)
    plotr( ppeak_struct.max_iteration.shifts(i,:),ppeak_struct.max_iteration.ints(i,:),'.','Color', [206, 194, 136]./255 ,'MarkerSize',9)
    hold on
end

highlightLine(metadata)

legend off

text( sort(mean(bucket_list,2)), ones(1, length(bucket_list)).*(abs(max(matrix(:)))./300.*-1), string(ppeak_struct.align_scale.metric),'HorizontalAlignment', 'left','Rotation',270)

hold on

legend( 'spectra' )
 h=plotr ( [nan,nan], '.','Color',  [206, 194, 136]./255,'MarkerSize',24,'DisplayName', char('Max of picked peaks'));
 h1=plotr ( [nan,nan], '-','Color', [111, 134, 149]./255,'LineWidth',2,'DisplayName', char('Bucket Boundaries'));
 h2=plotr ( [nan,nan], 'o','Color', [242, 112, 89]./255,'DisplayName', char('All picked peaks'));
 h3=plotr ( [nan,nan], '.','Color', [51, 108, 215]./255,'MarkerSize',24, 'DisplayName',char('Gap filled and max picked peaks'));



 title({'Alignment scale per bin and picked peaks'})
 
%%


%% Plotting the results in summary form 
% An overall view of the alignment scale, with buckets on the x axis and
% alignment scale on the yaxis
figure,
hold
plotr(1:length(ppeak_struct.align_scale.metric), (ppeak_struct.align_scale.metric));
%plot(1:length(ppeak_struct.align_scale.metric), sort(ppeak_struct.align_scale.metric), 'LineWidth',2)

title( {'Alignment metric summary'})
ylabel({'Alignment metric '}, 'FontSize',14)
xlabel({'Bucket Number'}, 'FontSize',15)
legend ([{'Metric per bucket'}, {'Sorted metric per bucket'}])
xticks(1:length(ppeak_struct.align_scale.metric))
xticklabels ( string(ppeak_struct.max_iteration.gap_f_median_ppm))
xtickangle(90)
%

%% Bar chart of the number of elements in each of the ranges of the alignment scale 

metric = ppeak_struct.align_scale.metric;

b_pl_5 = numel(metric (metric <= 5))./457*100;
b1 = ones(numel(metric (metric <= 5)),1);

b_pl_10 = numel(metric (metric >5 & metric<= 10  ))./457*100;
b2 =ones(numel(metric (metric >5 & metric<= 10  )),1).*2;

b_pl_15 = numel(metric (metric >10 & metric <= 15))./457*100;
b3 =ones(numel(metric (metric >10 & metric<= 15  )),1).*3;

b_pl_20 =numel( metric (metric >15 & metric <= 20))./457*100;
b4 =ones(numel(metric (metric >15 & metric<= 20  )),1).*4;

b_pl_25 = numel(metric (metric >20 & metric <= 25  ))./457*100;
b5 =ones(numel(metric (metric >20 & metric<= 25  )),1).*5;

b_pl_30 = numel(metric (metric >25 & metric <= 30  ))./457*100;
b6 =ones(numel(metric (metric >25 & metric<= 30  )),1).*6;

b_pl_40 = numel(metric( metric >30 & metric <= 40  ))./457*100;
b7 =ones(numel(metric (metric >30 & metric<= 40  )),1).*7;

b_pl_40plus = numel(metric( metric >40))./457*100;
b8 =ones(numel(metric (metric >40 )),1).*8;

bpl = [ b_pl_5; b_pl_10; b_pl_15; b_pl_20; b_pl_25; b_pl_30; b_pl_40; b_pl_40plus];


cc = [219, 254, 135; ...
    255, 227, 129;...
    206, 194, 136; ...
    111, 134, 149;...
    28, 68, 142; ...
    242, 112, 89;...
    83, 19, 30;...
    90, 70, 76]./255;

figure,
hold


b= bar (bpl,'FaceColor','flat');
for k = 1:size(bpl,1)
    b.CData(k,:) = cc(k,:);
end
xticks(1:8)
xticklabels ([ {'<5 - Excelent'}, {'5-10 - V. good'}, {'10-15 - Good'}, {'15-20 - OK'}, {'20-25 - Caution'}, {'25-30 - Inspect'}, {'30-40 - Not great'},{'40> - Unaligned'}])
xtickangle(15)
ylabel({'Percent of total number of peaks'}, 'FontSize',14)
xlabel({'Alignment metric ranges'}, 'FontSize',15)
title({'Defining a quality of alignment to scale values ranges'})












