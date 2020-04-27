function [VIP_Table,VIP_Sig_Table,VIP_Peaks_Sorted] = VIP_table(VIP_IN,Component_Number,Peaks_table,ppm,VIP_cutoff, intensity_cutoff)
%
%   Generates Three tables which present the information from VIP analysis of the PLS data. 
%       VIP_Table 
%            Contains VIP scores for all ppm values.
%       
%       VIP_Peaks_Sorted 
%           Truncates the VIP_Table to only show VIP scores corresponding to spectral peaks in the dataset.
%       
%       VIP_SIG_Table 
%           Lists only the VIP scores corresponding to ppm values from VIP_Peak_Sorted, where the VIP score satisfies a certain cutoff (VIP_cutoff)
%
%   Arguments:
%
%       VIP_IN
%           VIP score input for you PLS model; generated as the output from
%           using the 'vip' function from PLS Toolbox (Eigenvector)
%
%       Componet_Number
%           The PLS component of interest for generating a table of VIP
%           scores
%
%       Peaks_table
%           The normalized spectral table used to generate the PLS model
%
%       ppm
%           The list of ppm values corresponding to the intensity values
%           listed in the normalized spectral table used to generate the
%           PLS model
%
%       VIP_cutoff
%           The desired VIP score cutoff used to generate the truncated VIP_Sig_Table
%
%       intensity_cutoff
%           The spectral intensity cutoff, below which a spectral peak is
%           not considered a real spectral peak (this is essentially a
%           cutoff below which the spectral data is considered baseline and
%           not real measurements of spectral peaks)
%
VIP_Table = [ppm; transpose(VIP_IN(:,Component_Number));mean(Peaks_table)];

marker = 1;
for i = 2:length(VIP_Table)-1
    if VIP_Table (3,i) > intensity_cutoff && VIP_Table (3,i) > VIP_Table (3,i-1) && VIP_Table (3,i) > VIP_Table (3,i+1)
        VIP_Peaks(1,marker) = VIP_Table (1,i);
        VIP_Peaks(2,marker) = VIP_Table (2,i);
        VIP_Peaks(3,marker) = VIP_Table (3,i);
        marker = marker +1;
    end
end
%Only keep the rows that have a VIP equal to or greater than the set value
%above
clear i;
marker=1;
for i = 1:length(VIP_Peaks)
      if VIP_Peaks(2,i)>=VIP_cutoff
         VIP_Sig_Table(1,marker) = VIP_Peaks(1,i);
         VIP_Sig_Table(2,marker) = VIP_Peaks(2,i);
         VIP_Sig_Table(3,marker) = VIP_Peaks(3,i);
         marker = marker+1;
      end
end

VIP_Peaks_Sorted = flipud(sortrows(transpose (VIP_Peaks), 2));

clear i;
figure,
plotr(ppm, mean(Peaks_table,1));
hold on;
for i = 1:length(VIP_Sig_Table)
     plot(VIP_Sig_Table(1,i),VIP_Sig_Table(3,i),'r*');
end
hold off;