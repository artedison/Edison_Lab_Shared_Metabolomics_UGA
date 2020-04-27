function [S_Table]= STOCSY_Table (X,ppm,peaks_table,intensity_cutoff,corr_cutoff)

% Arguments    ppm = string of ppm values used as input for STOCSY Program
%   
%              X = the X value input used for the STOCSY Program
%
%              intensity_cutoff = the cutoff used to distinguish between
%                                 background intensity and real peaks
%
%              peaks_table = the table of significant peaks to use when
%              calculating STOCSY correlations for peaks
%
%              corr_cutoff = correlation value cutoff when finding peaks
%              that correlate with peak input used in STOCSY program
%
%
% Output       Structure table with the following cells:
%                An identifier called 'peakXXX' denoting the peak used for
%                STOCSY analysis, containing the table with the following
%                   Column 1. ppm value of correlated peak
%                   Column 2. average intensity of the correlated peak
%                   Column 3. correlation value of the peak
%                   Column 4. covariance value for the peak
%
%
%                Note, this table is sorted from highest corelation value peak to the lowest  

for i=1:length(peaks_table)
[A, B]=STOCSY_AJB(peaks_table(1,i),X,ppm);
[C] = Sig_STOCSY(ppm,X,intensity_cutoff,A,B,corr_cutoff);
ident = num2str(peaks_table(1,i));
ident = strrep(ident,'.','_');
ident = strrep(ident,'-','neg');
p = 'peak';
pident = strcat(p,'_',ident);
S_Table.(pident(1:end))= flipud(sortrows(transpose (C), 3));
end
