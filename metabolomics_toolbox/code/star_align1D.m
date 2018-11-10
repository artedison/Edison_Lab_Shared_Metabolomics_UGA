function XAL=star_align1D(X,ppm,represent,alignment_method,Seg_ppm,MaxShift_ppm,slack_ppm)

% XAL=star_align1D(X,ppm,represent,alignment_method,Seg_ppm,MaxShift_ppm,slack_ppm)
% 
% Calculates star alignment to a representative spectrum using
% alignment method of your choice - default is PAFFT
%
% Arguments:
% 
% X                      Data matrix of spectra
% ppm                    chemical shift vector
% represent              representative spectrum of spectral set- can be
%                        'mean','median','max', 'min', 'var', or an interger for the index of
%                        the spectrum in the full stack
% alignment_method       string, either 'CCOW','PARCCOW','ICOSHIFT','RAFFT','PAFFT'
%
% Optional Arguments:
%
% Seg_ppm                Length of segment in ppm (dafault 0.08)
% MaxShift_ppm           Maximum distance a segment can be shifted in ppm
%                        (default is 0.05)
% slack_ppm              Slack parameter for CCOW in ppm- default is 0.005 
%
% Return Values:
% XAL                    Star-aligned spectra
%
% Steven L Robinette


if isstr(represent)==0
    reference=X(represent,:);
else
    represent=str2func(represent);
    reference=represent(X);
end

if exist('slack_ppm')~=1
    slack=round(0.005/(ppm(2)-ppm(1)));
else
    slack=round(slack_ppm/(ppm(2)-ppm(1)));
end
if exist('MaxShift_ppm')~=1
    MaxShift=round(0.05/(ppm(2)-ppm(1)));
else
    MaxShift=round(MaxShift_ppm/(ppm(2)-ppm(1)));
end
if exist('Seg_ppm')~=1
    SegLength=round(0.08/(ppm(2)-ppm(1)));
else
    SegLength=round(Seg_ppm/(ppm(2)-ppm(1)));
end

NumSegs=round(size(X,2)/SegLength);
    
switch upper(alignment_method);
    case('CCOW')
        XAL=CCOW(X,reference,'SegLength',SegLength,'maxPeakShift',MaxShift,'Slack',slack);
        case('PARCCOW')
        XAL=parCCOW(X,reference,'SegLength',SegLength,'maxPeakShift',MaxShift,'Slack',slack);
    case('ICOSHIFT')
        XAL=icoshift(reference,X,NumSegs,MaxShift);
    case('RAFFT')
        XAL=RAFFT(X, reference, MaxShift, 1);
    case('PAFFT')
        XAL=PAFFT(X, reference, SegLength, MaxShift);
end

end