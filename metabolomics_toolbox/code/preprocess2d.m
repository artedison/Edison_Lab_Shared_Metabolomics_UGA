path{1}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/N2-1-cosy.fid/1';
path{2}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/N2-3-cosy.fid/1';
path{3}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/N2-5-cosy.fid/1';
path{4}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/N2-7-cosy.fid/1';
path{5}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/m130-9-cosy.fid/1';
path{6}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/m130-11-cosy.fid/1';
path{7}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/m130-13-cosy.fid/1';
path{8}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/m130-15-cosy.fid/1';
path{9}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/ok693-17-cosy.fid/1';
path{10}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/ok693-19-cosy.fid/1';
path{11}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/ok693-21-cosy.fid/1';
path{12}='/Users/sr2408/Desktop/Steve/Grad/Edison_Schroeder/Yev_Magn_spectra/ok693-23-cosy.fid/1';

cd('/Users/sr2408/Desktop/Steve/Current_Code/MetaboToolbox/2D')

Y=[1 1 1 1 0 0 0 0 0 0 0 0]';
Y2=[1 1 1 1 2 2 2 2 3 3 3 3]';

%Load Spectra
spectra=Load2D(path');

%Make 3D Matrix of 2D Spectra
Setup2D(spectra);

%check alignment
Y=rand(size(X,3));
stackplot(X,XNoise,ppm1,ppm1,Y)

%Segment spectra
label=segment2D(X,XNoise,ppm1,ppm2);

%Align if necessary
XAL=star_align2D(X,label,ppm1,ppm2,maxshift);

%Numerically integrate segments
binmat=bin2D(X,label);

normcheck(binmat)
normbin=normalize(binmat,'PQN');

varcheck(normbin)

XSN=scale(normbin,'logoff');

[sample_order,variable_order]=two_way_cluster(XSN,'weighted','spearman',Y);

[p,sigs]=MWAS(XSN,Y,'bonferroni');
manhattan(XSN,Y,1:size(XSN,2),p,sigs)

PCA=nipalsPCA(XSN,5);
VisScores(XSN,PCA,[1 2 3],Y2);

PLS=pls2D(Y,XSN,5,'da',10);
VisScores(XSN,PLS,[1 2],Y2);



