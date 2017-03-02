function alignedSpectrum = RAFFT(spectra, reference, shift, lookahead)

%	FUNCTION: alignedSpectrum = RAFFT(spectra, reference, shift, lookahead)
%
%	Spectrum alignment of spectral data to a reference using fast fourier
%	transform correlation theorem by recursive alignment
%	
%	INPUT:	spectra	    -	Spectra to be shift corrected (D x N - where
%                           D is the number of samples, and N is number of
%                           data points.
%			reference   -	Reference spectra, (spectrum and target must be
%			                of the same length).
%           shift       -   Optional: limits maximum shift allowed for each
%                           segment.
%           lookahead   -   Optional: Allows the recursive algorithm to
%                           look-a-head to check to local misalignments.
%   OUTPUT:	
%	alignedSpectrum - resulting aligned spectra
%
%MAIN----------------------------------------------------------------------
if length(reference)~=length(spectra)
    error('Reference and spectra of unequal lengths');
elseif length(reference)== 1
    error('Reference cannot be of length 1');
end
if nargin==2
    shift = length(reference);
    lookahead=1;
elseif nargin==3
    lookahead=1;
elseif nargin==4
    if lookahead==0
        lookahead=1;
    end
end

for i=1:size(spectra,1)
    aligned=recurAlign(spectra(i,:),reference,shift,lookahead);
    alignedSpectrum(i,:)= aligned;
end

%Recursive segmentation----------------------------------------------------
function aligned = recurAlign(spectrum, reference, shift, lookahead)
%stop when length of segment is less than 10
if length(spectrum) < 10
    aligned = spectrum;
    return;
end

lag = FFTcorr(spectrum,reference,shift);
%stop if the segment is perfectly aligned and there is no need to lookahead
if (lag == 0 && lookahead <= 0)
    aligned = spectrum;
    return;
else
    if lag == 0
        lookahead = lookahead-1;
    end
    aligned = spectrum;
    if abs(lag) < length(spectrum)
        aligned = move(spectrum,lag);
    end
    mid = findMid(aligned);
    firstSH= aligned(1,1:mid);
    firstRH= reference(1,1:mid);
    secSH = aligned(1,mid+1:length(aligned));
    secRH = reference(1,mid+1:length(reference));
    aligned1 = recurAlign(firstSH,firstRH,shift, lookahead);
    aligned2 = recurAlign(secSH,secRH,shift, lookahead);
    aligned = [aligned1 aligned2];
end

% FFT cross-correlation------------------------------------------------
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

%if the max segment correlation is very poor then assume no correlation and
%no lag.
if maxi < 0.1
    lag =0;
    return;
end
if maxpos > length(vals)/2
   lag = maxpos-length(vals)-1;
else
   lag =maxpos-1;
end

%shift segments----------------------------------------------------------
function movedSeg = move(seg, lag)
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
function mid = findMid(spec)

M=ceil(length(spec)/2);
specM=spec(1,M-floor(M/4):M+floor(M/4));
[C,I]=min(specM);
mid = I(1,1)+M-floor(M/4);