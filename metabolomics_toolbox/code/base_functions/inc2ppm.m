function peak_list_ppm = inc2ppm(header_data)

% ppm=inc2ppm(header_data);
% Greg Stupp

% Modified from:
% Covariance NMR Toolbox, ver. 1.0b, (C) (2010) David A. Snyder(1) along with Timothy Short(1),
%     Leigh Alzapiedi(1) and Rafael Br�schweiler (2)
%     (1) Department of Chemistry, College of Science and Health, William Paterson University
%     (2) Department of Chemistry and NHMFL, Florida State University
%
% License
%
% Copyright (c) 2010, David Snyder
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% * Redistributions of source code must retain the above copyright notice, this
%  list of conditions and the following disclaimer.
%
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% Currently, Cite As
% David Snyder (2022). Covariance NMR Toolbox
% (https://www.mathworks.com/matlabcentral/fileexchange/27264-covariance-nmr-toolbox)
% MATLAB Central File Exchange. Retrieved February 11, 2022.


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
