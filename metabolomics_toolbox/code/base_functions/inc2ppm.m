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
