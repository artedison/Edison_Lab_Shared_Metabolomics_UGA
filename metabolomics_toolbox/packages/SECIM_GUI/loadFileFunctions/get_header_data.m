function header_data = get_header_data(header_in)
% header_data = get_header_data(header_in)
%    reorganize key data in NMR pipe header into a matrix with each row having
%    a different data type and each column associated with a specific dim

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
