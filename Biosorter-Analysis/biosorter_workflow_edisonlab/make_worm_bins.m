function make_worm_bins(lnExt)
%% setting bin maximum values - logarithmic units
bin1 = -1.4;
bin2 = 0.4;
bin3 = 1.4;
bin4 = 2.3;
bin5 = 4;

%% break the maxtrix of log normalized extinction values into bins

for j = 1:length(lnExt);
TP_1_idx{1,j} = find(lnExt{1,j}(:,1) >= bin1 & lnExt{1,j}(:,1) < bin2);
TP_2_idx{1,j} = find(lnExt{1,j}(:,1) >= bin2 & lnExt{1,j}(:,1) < bin3);
TP_3_idx{1,j} = find(lnExt{1,j}(:,1) >= bin3 & lnExt{1,j}(:,1) < bin4);
TP_4_idx{1,j} = find(lnExt{1,j}(:,1) >= bin4 & lnExt{1,j}(:,1) < bin5);
end
%% count how many items in each bin
for w = 1:length(lnExt);
      TP_1_length1{1,w} = length(TP_1_idx{1,w});
      TP_2_length1{1,w} = length(TP_2_idx{1,w});
      TP_3_length1{1,w} = length(TP_3_idx{1,w});
      TP_4_length1{1,w} = length(TP_4_idx{1,w});
end
%% change to matrix and invert dimension
TP_1_length1 = cell2mat(TP_1_length1);
TP_2_length1 = cell2mat(TP_2_length1);
TP_3_length1 = cell2mat(TP_3_length1);
TP_4_length1 = cell2mat(TP_4_length1);

TP_1_length1 = TP_1_length1';
TP_2_length1 = TP_2_length1';
TP_3_length1 = TP_3_length1';
TP_4_length1 = TP_4_length1';

%% This matrix "m" gives the number of events in each bin. This is what is
% finally used to generate the histogram
m = [TP_1_length1 TP_2_length1 TP_3_length1 TP_4_length1];

%% plot histograms
prompt1 = 'How many rows in histogram subplot?';
prompt2 = 'How many columns in histogram subplot?';

subplotrows = input(prompt1,'s');
subplotcol = input(prompt2,'s');
%%
subplotrows = str2num(subplotrows);
subplotcol = str2num(subplotcol);
%%
[rows, ~] = size(m);
for k = 1:rows;
    subplot(subplotrows,subplotcol,k)
    bar(m(k,:));
end
end