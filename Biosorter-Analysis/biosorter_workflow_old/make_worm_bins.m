%% setting bin maximum values - logarithmic units
bin2 = 0.4;
bin3 = 1.4;
bin4 = 2.3;
bin1 = -1.4;
bin5 = 4;

%% break the maxtrix of normalized extinction values into bins

Ext4 = lnExt;

for j = 1:length(Ext4);
TP_1_idx{1,j} = find(Ext4{1,j}(:,1) >= bin1 & Ext4{1,j}(:,1) < bin2);
TP_2_idx{1,j} = find(Ext4{1,j}(:,1) >= bin2 & Ext4{1,j}(:,1) < bin3);
TP_3_idx{1,j} = find(Ext4{1,j}(:,1) >= bin3 & Ext4{1,j}(:,1) < bin4);
TP_4_idx{1,j} = find(Ext4{1,j}(:,1) >= bin4 & Ext4{1,j}(:,1) < bin5);
end

for w = 1:length(Ext4);
      TP_1_length1{1,w} = length(TP_1_idx{1,w});
      TP_2_length1{1,w} = length(TP_2_idx{1,w});
      TP_3_length1{1,w} = length(TP_3_idx{1,w});
      TP_4_length1{1,w} = length(TP_4_idx{1,w});
end

TP_1_length1 = cell2mat(TP_1_length1);
TP_2_length1 = cell2mat(TP_2_length1);
TP_3_length1 = cell2mat(TP_3_length1);
TP_4_length1 = cell2mat(TP_4_length1);

TP_1_length1 = TP_1_length1';
TP_2_length1 = TP_2_length1';
TP_3_length1 = TP_3_length1';
TP_4_length1 = TP_4_length1';

%% This matrix gives the number of events in each bin. This is what is
% finally used to generate the histogram
m = [TP_1_length1 TP_2_length1 TP_3_length1 TP_4_length1];
