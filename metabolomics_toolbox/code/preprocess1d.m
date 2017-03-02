%Load Spectra
path='/Users/sr2408/Desktop/Steve/Grad/Gahl/slr_MS_0311';
spectra=Load1D(path,'bruker');
spectra(158:end)=[];

%make 2D matrix of 1D spectra
[X,ppm,XTitles]=Setup1D(spectra);
figure, plot(ppm,X)

%remove unnecessary regions, check alignment
XR=remove_region(X,ppm,4.5,5);
figure, plot(ppm,XR)

%alignment
XALg=guided_alignment(XR,ppm);
XAL=star_alignment(XR,ppm,'mean','PAFFT');

%peak picking
[peakmatrix,shifts]=Peakpick1D(XALg,ppm,'mean',.9);