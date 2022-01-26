function [ppeak_struct] = peaks_per_bin (matrix, ppm,buckets, ppeak_struct);

% Requirements:
%   This function requires the use of Peakpick1D_per_spectrum structure output
%   Requires input of buckets from optBucket function

% Summary:
% This operates on the buckets and removes bins with no width (boundaries
% are the same) and sorts them from smallest to biggest ppm
% 
% This function is an operation matrix function, it simply uses the boundaries defined by
% bucket_list to find:
% the number of picked peaks from Peakpick1D_per_spectra output,
% it then calculates the maximum intensity peak within each bucket and
% its respective ppm shift (max_iteration)
% while also counting the number of peaks in each of the max iterations across
% all the spectra.
% finally it finds the ppms and intensities where there were no peaks
% detected per spectrum and then assings the minimum value of the raw
% data and the median ppm for that particular bucket (gap filling)


% input:
%   matrix  - -  where rows are samples and columns are ppms
%   ppm      - -  ppm vector with the same number of columns as matrix and
%                                    one row
%   buckets - - A matrix of 2 columns where each row specifies one bucket
%   ppeak_struct - - A structure from Peakpick1D_per_spectra containing a
%   field 'iteration' and inside 'iteration' contains 2 fields 'ints' and
%   'shifts'
%
%
% Output:
%
%   updates the structure ppeak_struct with the following fields:
%   ppeak_struct.max_iteration - containing -
%   ppeak_struct.max_iteration.sorted_bucket_list - This is the bucket list
%                                                                             sorted and boundaries that do not contain peaks removed
%                                                                  ppeak_struct.max_iteration.ints - Intensity per spectra per bucket
%                                                                  ppeak_struct.max_iteration.shift - The ppm value per bucket per spectra
%                                                                  ppeak_struct.max_iteration.max_elements_per_bin - number of peaks per bin over all the spectra
%                                                                  ppeak_struct.max_iteration.gap_f_shift
%                                                                  ppeak_struct.max_iteration.gap_f_ints
%                                                                  ppeak_struct.max_iteration.gap_f_median_ppm
%
%    ppeak_struct.prominence - containing - min_shift - The minimum matrix value
%                                                                                     per bin
%    ppeak_struct.prominence - containing -  min_ints - The minimum matrix value
%                                                                                     ppm value per bin
%
%    ppeak_struct.prominence.prominence_ints - the max peak picked gap filled
%                                                                       intensity minus the minimum matrix value for that bin and that spectra
% 
% 
% 

% Note: Structures names can be user defined, however substructures and
% fields need to remain as defined by each function 

%%
%% Sorting the refine bucket list so its from bigger to smaller,

ppeak_struct.max_iteration.sorted_bucket_list = sort( buckets );

bucket_list = ppeak_struct.max_iteration.sorted_bucket_list;

%% ensures the first element is the smallest bin boundary and the second element is the largest 

for i = 1:length(ppeak_struct.max_iteration.sorted_bucket_list )

if bucket_list(i,1)>bucket_list(i,2) == 1
    bucket_list(i,1)= ppeak_struct.max_iteration.sorted_bucket_list(i,2);
    bucket_list(i,2)= ppeak_struct.max_iteration.sorted_bucket_list(i,1);
    
else
    bucket_list(i,1) = bucket_list(i,1);
    bucket_list(i,2) = bucket_list(i,2);
    
end
end

bucket_list(bucket_list (:,1)==bucket_list(:,2),:) = []; % removes bins that have the same value for both starting and finnishing boundaries

%% Calculating the number of peaks in each spectra given the bins defined by the optbucket


for i = 1:size(matrix,1)
    for    ii=1:length (bucket_list )
        
        ppeak_struct.iteration(i).total_elements_per_bin(ii) = numel(...
            ppeak_struct.iteration(i).shifts(:,...
            ppeak_struct.iteration(i).shifts>bucket_list (ii,1)...
            &  ppeak_struct.iteration(i).shifts<bucket_list (ii,2) ...
            ) ...
            );
        
      echo off  
    end
end


%% Calculating the max peak from the picked peak data

for i = 1:size(matrix,1)
    for    ii=1:length (bucket_list)
        
        if  isempty(ppeak_struct.iteration(i).ints(:,...
                ppeak_struct.iteration(i).shifts>=bucket_list(ii,1)...
                &  ppeak_struct.iteration(i).shifts<=bucket_list(ii,2) ) )  == 1 ; %if there are no peaks between the bondaries add 'nan' as a value
        
            ppeak_struct.max_iteration.ints(i,ii)= 0;
            ppeak_struct.max_iteration.shifts(i,ii)= 0;
            
        else
            [ ppeak_struct.max_iteration.ints(i,ii),pp(i,ii)] = ... 
                max( ppeak_struct.iteration(i).ints(:,...
                ppeak_struct.iteration(i).shifts>=bucket_list(ii,1)...
                &  ppeak_struct.iteration(i).shifts<=bucket_list(ii,2) )...
                ,[],2) ;
            
            ppmp = ppeak_struct.iteration(i).shifts (...
                ppeak_struct.iteration(i).shifts>=bucket_list(ii,1)...
                &  ...
                ppeak_struct.iteration(i).shifts<=bucket_list(ii,2) );
            
            ppeak_struct.max_iteration.shifts(i,ii) = ppmp(pp(i,ii));
            
            
            

        end
        
    end
  echo off  
end
clearvars ppmp pp

%% Clearing columns that have no peaks in both ppm, bucket list and int matrix

% ppeak_struct.max_iteration.shift = ppeak_struct.max_iteration.shift(:,~all(isnan(ppeak_struct.max_iteration.ints)));
% 
% ppeak_struct.max_iteration.sorted_bucket_list = ppeak_struct.max_iteration.sorted_bucket_list (~all(isnan(ppeak_struct.max_iteration.ints))',:);
% 
% ppeak_struct.max_iteration.ints = ppeak_struct.max_iteration.ints(:,~all(isnan(ppeak_struct.max_iteration.ints)));
% 
% 
% bucket_list = ppeak_struct.max_iteration.sorted_bucket_list ;
%% calculating the number of peaks per bin in the max peak picked matrix


% for    i=1:length (bucket_list )
%     
%     ppeak_struct.max_iteration.max_elements_per_bin(:,i) = numel(...
%         ppeak_struct.max_iteration.shift(ppeak_struct.max_iteration.shift(:,i)>0, i));
%   
%     echo off
% end

%% Gap filling ppms for bins that had peaks that were not detected

for i= 1:size(bucket_list,1)
    
    zero_peaks =ppeak_struct.max_iteration.shifts(:,i)==0; %find the ppms of peaks that were not picked by Peakpick1D_per_spectra
    
    if  sum(int8(zero_peaks))<= size(matrix,1) % checks for empty columns
        
        ppeak_struct.max_iteration.gap_f_shifts(:,i)=ppeak_struct.max_iteration.shifts(:,i);
        
        ppeak_struct.max_iteration.gap_f_shifts(...
            zero_peaks,...
            i) = median(ppeak_struct.max_iteration.shifts(ppeak_struct.max_iteration.shifts(:,i)>0 | ppeak_struct.max_iteration.shifts(:,i)<0,i));
        
    else
        ppeak_struct.max_iteration.gap_f_shifts(:,i) = ppeak_struct.max_iteration.shifts(:,i);
        disp(['empty column',string(i)])
        
    end
    echo off
end



%%  %% Gap filling intensities for bins that had peaks that were not detected

for i= 1:size(bucket_list,1)
    
    nan_peaks =isnan(ppeak_struct.max_iteration.ints(:,i));
    
    if  sum(int8(nan_peaks))<= size(matrix,1)
        
        ppeak_struct.max_iteration.gap_f_ints(:,i)=ppeak_struct.max_iteration.ints(:,i);
        
        ppeak_struct.max_iteration.gap_f_ints(...
            isnan(ppeak_struct.max_iteration.ints(:,i)),...
            i) ...
            = min(matrix(:, matchPPMs(median(ppeak_struct.max_iteration.gap_f_shifts(:,i)),ppm) ));
        
    else
        ppeak_struct.max_iteration.gap_f_ints(:,i)=ppeak_struct.max_iteration.ints(:,i);
        disp(['empty column',string(i)])
    end
    echo off
end

ppeak_struct.max_iteration.sorted_bucket_list=bucket_list;

%% Getting a single ppm vector with median ppm for each bin/peak

 ppeak_struct.max_iteration.gap_f_median_ppm = median( ppeak_struct.max_iteration.sorted_bucket_list,2);
 
 %%  %% Finding the minimum intensity in the raw data for each bin per spectra and its shift 
% 
% for i= 1:size(matrix,1)
%     mat = matrix (i,:);
%     
%     for ii = 1:size(bucket_list,1)
%         
%        [ min_ints(:,ii), min_shift(:,ii)] =min(...
%            mat( :,...
%                 ppm>=bucket_list(ii,1)...
%                 &  ...
%                 ppm<=bucket_list(ii,2)...
%                 )...
%                 );
%             
%             pp = ppm( :,...
%                 ppm>=bucket_list(ii,1)...
%                 &  ...
%                 ppm<=bucket_list(ii,2)...
%                 );
%             
%             min_shift(:,ii) = pp(min_shift(:,ii) );
%     end
%      ppeak_struct.prominence.min_shift(i,:) = min_shift;
%      ppeak_struct.prominence.min_ints(i,:) = min_ints;
%      
% end


 %%   Calculating the max peak picked prominence for each spectra 
%  for i = 1:size(bucket_list,1)
%      
%      ppeak_struct.prominence.prominence_ints (:,i) = ppeak_struct.max_iteration.gap_f_ints(:,i) -  ppeak_struct.prominence.min_ints(:,i);
% 
%   ppeak_struct.prominence.prominence_ints(ppeak_struct.prominence.prominence_ints (:,i) <=0,i) = ppeak_struct.max_iteration.gap_f_ints (ppeak_struct.prominence.prominence_ints (:,i) <=0,i);
%  end






