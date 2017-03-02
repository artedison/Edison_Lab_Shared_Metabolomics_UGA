function spectra=LoadPipe(name)
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

if (nargin==0)
    name=strtrim(ls('*.ft'));
    title=name(1:length(name)-3);
else
    [pathstr, name2, ext] = fileparts(name);
    if length(pathstr)>1
        title=name2;
    else
        title=name(1:length(name)-3);
    end
end


[spectra.real axes]=read_nmrp(name);

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


function [spec_array axes_info] = read_nmrp(input_fname,donor_dims)
% [spec_array axes_info] = read_nmrp(input_fname)
% [spec_array axes_info] = read_nmrp(input_fname,donor_dims)
%     read in NMR pipe format file in input fname, parse the header and
%     return an array as well as a matrix tabulating what info is in which
%     axis.  Note that the dimension order of NMR pipe format spectra is
%     w4 w3 w2 w1 = dims # 1, 2, 3, 4.  The donor dims (i.e. the dims to be
%     subsumed by covariance) can be optionally input as a vector, e.g.
%     [1 3 4].  If the donor dims are not indicated, the first N/2 dims of 
%     the input N-D spectrum are considered acceptor dims and the remaining 
%     dims are considered donor dims to be covaried.

% Covariance NMR Toolbox, ver. 1.0b, (C) (2010) David A. Snyder(1) along with Timothy Short(1),
%     Leigh Alzapiedi(1) and Rafael Br�schweiler (2)
%     (1) Department of Chemistry, College of Science and Health, William Paterson University
%     (2) Department of Chemistry and NHMFL, Florida State University
    

number_args_in = nargin;
    
spectrum_in = fopen(input_fname, 'r');

header_in = fread(spectrum_in, 512, 'float32');
axes_info = get_header_data(header_in);
size_row = max(axes_info(26,:),1); % sometimes unused dims may have size set to 0, we need these set to 1
size_spec = prod(size_row); % calculate size of spec 

data = fread(spectrum_in, size_spec, 'float32');
spec_array = squeeze(reshape(data,size_row));

clear data;

% add a row to axes_info indicating donor/acceptor dims

num_dims = sum(axes_info(26,:) > 1);

if(number_args_in == 2) % we have a given donor dim set, so consider acceptors as anything else
	axes_info(27,:) = -1; % -1 is code for unused dim (which is == to donor dim)
    axes_info(27,1:num_dims) = 1;
	axes_info(27,donor_dims) = -1; % -1 is code for donor dims	
elseif(number_args_in == 1)
	num_dims = sum(size_row > 1);
	halfnum_dims = floor(num_dims/2);

	axes_info(27,:) = [ones(1,halfnum_dims) -1*ones(1,4-halfnum_dims)];
end
end

function peak_list_ppm = inc2ppm(header_data)

% ppm=inc2ppm(header_data);
% Greg Stupp

% Modified from:
% Covariance NMR Toolbox, ver. 1.0b, (C) (2010) David A. Snyder(1) along with Timothy Short(1),
%     Leigh Alzapiedi(1) and Rafael Br�schweiler (2)
%     (1) Department of Chemistry, College of Science and Health, William Paterson University
%     (2) Department of Chemistry and NHMFL, Florida State University
    
    
% peak_list_ppm = inc2ppm(peak_list_inc, header_data)
%    map a peak list (in incriment units) to corresponding ppm values

Ni_row = header_data(26,:);

%added by GS
num_dims = sum(Ni_row > 1);

%from orig
%num_dims = min(size(peak_list_inc,2),sum(Ni_row > 1));
%num_peaks = size(peak_list_inc,1);

Ni_row = Ni_row(1:num_dims);
right_ppm = header_data(16,1:num_dims)./header_data(14,1:num_dims);
sw = header_data(20,1:num_dims)./header_data(14,1:num_dims);
left_ppm = right_ppm + sw;

%From original
%peak_list_ppm = (-repmat(sw,num_peaks,1)./repmat(Ni_row,num_peaks,1)).*peak_list_inc(:,1:num_dims) + repmat(left_ppm,num_peaks,1);

%Added by GS
for d=1:num_dims
    a=[right_ppm(d):(left_ppm(d)-right_ppm(d))/Ni_row(d):left_ppm(d)];
    a(end)=[];
    a=fliplr(a);
    a=a';
    eval(['peak_list_ppm.ppm' int2str(d) '=a;']);
end

end


function header_data = get_header_data(header_in)
% header_data = get_header_data(header_in)
%    reorganize key data in NMR pipe header into a matrix with each row having
%    a different data type and each column associated with a specific dim

% Covariance NMR Toolbox, ver. 1.0b, (C) (2010) David A. Snyder(1) along with Timothy Short(1),
%     Leigh Alzapiedi(1) and Rafael Br�schweiler (2)
%     (1) Department of Chemistry, College of Science and Health, William Paterson University
%     (2) Department of Chemistry and NHMFL, Florida State University

UNUSED_VALUE = -666; % FD's alt. code for zero: there seems to be no SIZE for y-dim nor LB for z and a-dims

header_data = [ header_in( 95 +1)  header_in( 428 +1)  header_in( 50 +1)  header_in( 53 +1) ; ... % APOD	1
   header_in( 413 +1)  header_in( 414 +1)  header_in( 400 +1)  header_in( 405 +1) ; ... % APODCODE		2
   header_in( 415 +1)  header_in( 420 +1)  header_in( 401 +1)  header_in( 406 +1) ; ... % APODQ1		3
   header_in( 416 +1)  header_in( 421 +1)  header_in( 402 +1)  header_in( 407 +1) ; ... % APODQ2		4
   header_in( 417 +1)  header_in( 422 +1)  header_in( 403 +1)  header_in( 408 +1) ; ... % APODQ3		5
   header_in( 64 +1)  header_in( 475 +1)  header_in( 476 +1)  header_in( 477 +1) ; ... % AQSIGN			6
   header_in( 418 +1)  header_in( 423 +1)  header_in( 404 +1)  header_in( 409 +1) ; ... % C1			7
   header_in( 66 +1)  header_in( 67 +1)  header_in( 68 +1)  header_in( 69 +1) ; ... % CAR				8
   header_in( 79 +1)  header_in( 80 +1)  header_in( 81 +1)  header_in( 82 +1) ; ... % CENTER			9
   header_in( 220 +1)  header_in( 222 +1)  header_in( 13 +1)  header_in( 31 +1) ; ... % FTFLAG			10
   header_in( 96 +1)  header_in( 98 +1)  header_in( 200 +1)  header_in( 201 +1) ; ... % FTSIZE			11
   header_in( 16 +1)  header_in( 18 +1)  header_in( 20 +1)  header_in( 22 +1) ; ... % LABEL			12
   header_in( 111 +1)  header_in( 243 +1)  UNUSED_VALUE    UNUSED_VALUE   ; ... % LB				13
   header_in( 119 +1)  header_in( 218 +1)  header_in( 10 +1)  header_in( 28 +1) ; ... % OBS			14
   header_in( 480 +1)  header_in( 481 +1)  header_in( 482 +1)  header_in( 483 +1) ; ... % OFFPPM		15
   header_in( 101 +1)  header_in( 249 +1)  header_in( 12 +1)  header_in( 30 +1) ; ... % ORIG			16
   header_in( 109 +1)  header_in( 245 +1)  header_in( 60 +1)  header_in( 62 +1) ; ... % P0			17
   header_in( 110 +1)  header_in( 246 +1)  header_in( 61 +1)  header_in( 63 +1) ; ... % P1			18
   header_in( 56 +1)  header_in( 55 +1)  header_in( 51 +1)  header_in( 54 +1) ; ... % QUADFLAG			19
   header_in( 100 +1)  header_in( 229 +1)  header_in( 11 +1)  header_in( 29 +1) ; ... % SW			20
   header_in( 386 +1)  header_in( 387 +1)  header_in( 388 +1)  header_in( 389 +1) ; ... % TDSIZE		21
   header_in( 152 +1)  header_in( 234 +1)  header_in( 58 +1)  header_in( 59 +1) ; ... % UNITS			22
   header_in( 257 +1)  header_in( 259 +1)  header_in( 261 +1)  header_in( 263 +1) ; ... % X1			23
   header_in( 258 +1)  header_in( 260 +1)  header_in( 262 +1)  header_in( 264 +1) ; ... % XN			24
   header_in( 108 +1)  header_in( 437 +1)  header_in( 438 +1)  header_in( 439 +1) ; ... % ZF			25
   header_in( 99 +1)  header_in( 219 +1)    header_in( 15 +1)  header_in( 32 +1) ]; % SIZE			26
end