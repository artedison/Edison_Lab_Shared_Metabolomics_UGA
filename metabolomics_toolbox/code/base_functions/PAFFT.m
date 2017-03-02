function alignedSpectrum = PAFFT(spectra, reference, segSize, shift)

%	FUNCTION: alignedSpectrum = PAFFT(spectra, target, segmentSize, shift)
%
%	Spectrum alignment of spectral data to a reference using fast fourier
%	transform correlation theorem by segmentation alignment
%	
%	INPUT:	spectra	    -	Spectra to be shift corrected (D x N - where
%                           D is the number of samples, and N is number of
%                           data points.
%			reference   -	Reference spectra, (spectrum and target must be
%			                of the same length).
%           segmentSize -   Minimum segement size allowed.
%           shift       -   Optional: limits maximum shift allowed for each
%                           segment.

%   OUTPUT:	
%	        alignedSpectrum - resulting aligned spectra
%MAIN----------------------------------------------------------------------
if length(reference)~=length(spectra)
    error('Reference and spectra of unequal lengths');
elseif length(reference)== 1
    error('Reference cannot be of length 1');
end
if nargin==3
    shift = length(reference);
end
for i=1:size(spectra,1)
    startpos = 1;
    aligned =[];
    while startpos <= length(spectra)
        endpos=startpos+(segSize*2);
        if endpos >=length(spectra)
            samseg= spectra(i,startpos:length(spectra));
            refseg= reference(1,startpos:length(spectra));
        else
            samseg = spectra(i,startpos+segSize:endpos-1);
            refseg= reference(1,startpos+segSize:endpos-1);
            minpos=findMin(samseg,refseg);
            endpos = startpos+minpos+segSize;
            samseg = spectra(i,startpos:endpos);
            refseg= reference(1,startpos:endpos);
        end
        lag = FFTcorr(samseg,refseg,shift);
        aligned = [aligned move(samseg,lag)];
        startpos=endpos+1;
    end
    alignedSpectrum(i,:)=aligned;
end

% FFT cross-correlation----------------------------------------------------
function lag = FFTcorr(spectrum, target, shift)
%padding
M=size(target,2);
diff = 1000000;
for i=1:20
    curdiff=((2^i)-M);
    if (curdiff > 0 && curdiff<diff)
        diff = curdiff;
    end
end

target(1,M+diff)=0;
spectrum(1,M+diff)=0;
M= M+diff;
X=fft(target);
Y=fft(spectrum);
R=X.*conj(Y);
R=R./(M);
rev=ifft(R);
vals=real(rev);
maxpos = 1;
maxi = -1;
if M<shift
    shift = M;
end

for i = 1:shift
    if (vals(1,i) > maxi)
        maxi = vals(1,i);
        maxpos = i;
    end
    if (vals(1,length(vals)-i+1) > maxi)
        maxi = vals(1,length(vals)-i+1);
        maxpos = length(vals)-i+1;
    end
end

if maxi < 0.1
    lag =0;
    return;
end
if maxpos > length(vals)/2
   lag = maxpos-length(vals)-1;
else
   lag =maxpos-1;
end
%-------------------------------------------------------------------

%shift segments----------------------------------------------------------
function movedSeg = move(seg, lag)
% sideways movement
if lag == 0 || lag >= length(seg)
    movedSeg = seg;
    return
end

if lag > 0
	ins = ones(1,lag)*seg(1);
	movedSeg = [ins seg(1:(length(seg) - lag))];
elseif lag < 0
	lag = abs(lag);
	ins = ones(1,lag)*seg(length(seg));
	movedSeg = [seg((lag+1):length(seg)) ins];
end

%find position to divide segment--------------------------------------
function minpos = findMin(samseg,refseg)

[Cs,Is]=sort(samseg);
[Cr,Ir]=sort(refseg);
minposA = [];
minInt = [];
for i=1:round(length(Cs)/20)
    for j=1:round(length(Cs)/20)
        if Ir(j)==Is(i);
            minpos = Is(i);
            return;
        end
    end
end
minpos = Is(1,1);