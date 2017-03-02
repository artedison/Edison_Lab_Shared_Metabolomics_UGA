function spectra=pipe2matlab(filepath)
% spectra=pipe2matlab(name)
%
% Greg Stupp UF
% Reads in an nmrPipe ft file. Outputs
% spectra.real,spectra.ppm1,spectra.ppm2,spectra.Title
% Uses functions borrowed from Covariance NMR Toolbox, ver. 1.0b, (C) (2010) David A. Snyder(1) along with Timothy Short(1),
%   Leigh Alzapiedi(1) and Rafael Br�schweiler (2)
%
% Arguments:
% name                 name of file (including extension)

[pathstr, title, ext] = fileparts(filepath);
[spectra.real axes]=read_nmrp(filepath);

ppm=inc2ppm(axes);
if length(fieldnames(ppm))==1
    spectra.ppm=ppm.ppm1;
else
    spectra.ppm1=ppm.ppm1;
    spectra.ppm2=ppm.ppm2;
end
spectra.real=spectra.real';
spectra.Title=title;
end


