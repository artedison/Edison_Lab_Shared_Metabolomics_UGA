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
    fclose(spectrum_in);    % MJ added 3JAN2018. Too many files left open without this; Matlab will crash with ~270 spectra
end
