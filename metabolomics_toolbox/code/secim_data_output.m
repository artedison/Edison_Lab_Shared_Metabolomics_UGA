%% output files setup

X_raw=vertcat(ppm,X_raw)';
X_removed=vertcat(ppm,X_removed)';
X_baseline_corrected=vertcat(ppm,X_baseline_corrected)';
X_aligned=vertcat(ppm,X_aligned)';
X_normalized=vertcat(ppm,X_normalized)';

%% write the output files

csvwrite('BR2_BLM_1_raw',X_raw);
csvwrite('BR2_BLM_1_removed',X_removed);
csvwrite('BR2_BLM_1_baseline_corrected',X_baseline_corrected);
csvwrite('BR2_BLM_1_aligned',X_aligned);
csvwrite('BR2_BLM_1_normalized',X_normalized);

%% Import the data 


[~, ~, BR2BLM1samplesheet] = xlsread('/Users/balm22/Desktop/SECIM_SB_NAMPT_blee/BR2_BLM_1/BR2_BLM_1_sample_info/BR2_BLM_1_sample_sheet.xlsx','Sheet1');
BR2BLM1samplesheet(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),BR2BLM1samplesheet)) = {''};


%% write the design_file

SampleID=BR2BLM1samplesheet(:,4);
group_info=BR2BLM1samplesheet(:,6);

design_file=horzcat(SampleID,group_info);

%% this is a new script.. you will need to pull from the git to get it

cell2csv('BR2_BLM_1_design_file2.csv',design_file);