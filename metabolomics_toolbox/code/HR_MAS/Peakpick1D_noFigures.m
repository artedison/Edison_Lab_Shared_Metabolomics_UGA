function [peakmatrix,shifts,p]=Peakpick1D_noFigures(X,ppm,represent,peakthresh,mode)

% [peakmatrix,shifts]=Peakpick1D(X,ppm,represent,peakthresh,mode)
% 
% Detects peaks in 1D NMR spectra and constructs matrix of peak intensities
% with chemical shifts of each peak maximum. 
%
% Arguments:
% X                    Data matrix of spectra
% ppm                  Chemical shift vector corresponding to X
% represent            representative spectrum of spectral set- can be
%                      'mean','median','max', 'min', 'var', or an interger for the index of
%                      the spectrum in the full stack
% peakthresh	       values between 0 and 1, higher threshold is greater
%                      specificity
% mode                 'Simple' (multiple of noise), or 'Complex'
%
% Return Values:
% peakmatrix           Matrix of peak intensities
% shifts               vector of chemical shifts at each position in
%                      peakmatrix

if exist('peakthresh')~=1
    peakthresh=.8;
end
if exist('mode')~=1
    mode='Complex';
end

p = reportParams('exclude',{'X','ppm'});

switch mode
    case 'Simple'
        if isstr(represent)==0
            spect=X(represent,:);
            pretargets=simplepick(spect,ppm,peakthresh);
        else
            represent=str2func(represent);
            spect=represent(X);
            pretargets=simplepick(spect,ppm,peakthresh);
        end
    case 'Complex'
        if isstr(represent)==0
            spect=X(represent,:);
            pretargets=PeakPicking(spect,ppm,peakthresh);
        else
            represent=str2func(represent);
            spect=represent(X);
            pretargets=PeakPicking(spect,ppm,peakthresh);
        end
end

% figure, plot(ppm,X)
% set(gca,'XDir','reverse')
% hold on
% for k=1:length(pretargets)
%     n(k)=scatter(ppm(pretargets(k)),max(X(:,pretargets(k))),'b');
% end

shifts=ppm(pretargets);
peakmatrix=X(:,pretargets);

end



function Validpeaks=PeakPicking(Sp,ppm,thresh)
% call configureRSPA to get peak peaking parameters

configureRSPA(ppm,thresh)
% perform Savitzkiy Golay smoothing
SpDerivs = sgolayDeriv(Sp,peakParam.iOrder,peakParam.iFrameLen,2);
SpSmooth = sgolayDeriv(Sp,peakParam.iOrder,peakParam.iFrameLen,1);
% indentify peaks
peaks=peakPeaks(SpSmooth,SpDerivs,Sp,[]);
% validate peaks
[peaksValidated,peakParam,minsegwidth]=validatePeaks(SpSmooth,peaks,peakParam,0);
Validpeaks=[peaksValidated.maxPos];
end

function configureRSPA(ppm,thresh)

%%%%%%%%%%%%%%%%%%%% Peak peaking parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
peakParam.thresh=thresh;
peakParam.ampThr = []; % amplitude value to threshold small peaks % 
peakParam.minPeakWidth = 0.005; %min peak width in ppm scale
peakParam.iFrameLen=0.005; %Savitzky-Golay frame length in ppm scale
peakParam.iOrder=3; %polynomial order of Savitzky - Golay filter
peakParam.peakEdgeMax=0.2; 


if ~isempty(ppm)
    peakParam.minPeakWidth=ppmToPt(peakParam.minPeakWidth,0,ppm(2)-ppm(1));
    peakParam.iFrameLen=ppmToPt(peakParam.iFrameLen,0,ppm(2)-ppm(1));
end

assignin('caller', 'peakParam',   peakParam);
return;
end

function  ampThr=getAmpThr(peaks,peakParam)
% Automatic determination of amplitude threshold for peak peaking
% based on the 5% of the most intensive peaks 

PeakCount=length(peaks);
peakMaxValues = repmat(NaN, [1, PeakCount]);

for i=1:PeakCount
      peakMaxValues(i)=peaks(i).maxVal-peaks(i).basl;
end

%%% SELECT THRESHOLD BASED ON (1-thesh)% OF THE MOST INTENSIVE PEAKS

index=floor(PeakCount*peakParam.thresh);
peakSortedValuess=sort(peakMaxValues);
ampThr=peakSortedValuess(index);
return;
end

function dpDerivs = sgolayDeriv(dpSpectr, iOrder,iFrameLen,j)
% Calculate smoothed derivates using Savitzky - Golay filter
% iFrameLen- the length of frame window

if nargin<1
    error('Incorrect number of input arguments');
end

if nargin<2 
    iOrder = 3; 
end

if nargin<3
    iFrameLen=11;
end

if nargin<4
    j=2; %Derivative
end

iFrameLen=(floor(iFrameLen./2))*2+1; % iFramLen must be odd

iSpecLen = length(dpSpectr);

[b,g] = sgolay(iOrder,iFrameLen);

dpDerivs(1:iFrameLen) = 0;
dpDerivs(iSpecLen-(iFrameLen+1)/2:iSpecLen) =0;

for n = (iFrameLen+1)/2:iSpecLen-(iFrameLen+1)/2
    %calculate first order derivate
    dpDerivs(n)=g(:,j)'*dpSpectr(n - (iFrameLen+1)/2 + 1: n + (iFrameLen+1)/2 - 1)';   
end

return;
end


function peaks=peakPeaks(SpSmooth,dpDerivs,Sp,debug)
% Peak peaking:
% Input: SpSmooth - smoothed spectrum
%        dpDerivs - smoothed derivative of the spectrum
%        debug - yes
% - the peak is identified if derivative crosses zero,
% i.e. sign(X'(i))>sing(X'(i+1))


if nargin<2
    % TODO: remove comment
    %error('Invalid number of input arguments')
end

iSpecLen=length(SpSmooth);
iPeakInd =1;

% Matrix pre-location
peaks(iSpecLen).maxPos=[];

for i=1:iSpecLen-1
    % coarse peak maximum position
    if dpDerivs(i)>=0&&dpDerivs(i+1)<0
        peaks(iPeakInd).maxPos=i+1;
        % Temporary starting and ending peak positions
        iPeakInd=iPeakInd+1;
    end
end

peakCount=iPeakInd-1;
peaks=peaks(1:peakCount);

targetPkIdx=1;

for srcPkIdx=1:peakCount
    maxPos=peaks(srcPkIdx).maxPos;
    
    while (maxPos > 2 && maxPos < iSpecLen-2)
        if SpSmooth(maxPos-1)<=SpSmooth(maxPos)&&...
                SpSmooth(maxPos)>=SpSmooth(maxPos+1)
            
           if(targetPkIdx > 1 && peaks(targetPkIdx-1).maxPos==maxPos)
                % the same maximum value - just skip it
                break;
            end
            % save the new index:
            peaks(targetPkIdx).maxPos = maxPos;
            targetPkIdx = targetPkIdx + 1;
            break;
        end
        if SpSmooth(maxPos)<=SpSmooth(maxPos+1)
            maxPos=maxPos+1;
        elseif SpSmooth(maxPos)<=SpSmooth(maxPos-1)
            maxPos=maxPos-1;
        end
    end
end

peakCount=targetPkIdx-1;
peaks=peaks(1:peakCount);

for i=1:peakCount
    j=peaks(i).maxPos;
    k=peaks(i).maxPos;

    % left boundary
    while SpSmooth(j)>=SpSmooth(j-1) && j-1~=1 %first index
        j=j-1;
    end

    % right boundary
    while SpSmooth(k)>=SpSmooth(k+1) && k+1~=iSpecLen %last index
        k=k+1;
    end

    peaks(i).startPos=j;
    peaks(i).endPos=k;
    peaks(i).centre=ceil((k+j)/2);
    peaks(i).startVal=SpSmooth(j);
    peaks(i).endVal=SpSmooth(k);
    peaks(i).index=i;
    % Use peak maximum position from original spectrum
    % instead of smoothed one.
    %peaks(i).maxVal=SpSmooth(peaks(i).maxPos);
    [peaks(i).maxVal, maxInd]=max(Sp(j:k));
    peaks(i).maxPos=j+maxInd-1;

    %estimate the baseline as minimum value:
    peaks(i).basl = min([SpSmooth(k), SpSmooth(j)]);
end

if ~isempty(debug);
    peakStartPositions=[];
    peakEndPositions=[];
    peakMaxPositions=[];
    for i=1:peakCount
        peakStartPositions=[peakStartPositions peaks(i).startPos];
        peakEndPositions=[peakEndPositions peaks(i).endPos];
        peakMaxPositions=[peakMaxPositions peaks(i).maxPos];
    end
    plotPeaksOrSegments(Sp,peakMaxPositions,peakStartPositions,...
        peakEndPositions,debug)
end

return;
end


function [validatedPeaks,peakParam,minsegwidth]=validatePeaks(SpSmooth,peaks,...
    peakParam,debug)
% input:          Peak peaking details
%                 peaks.  maxPos - peak maxium position
%                         startPos - start position
%                         endPos - end position
%                         maxVal - maximum value
%                         startVal - start value
%                         endVal - end value
%                         basl - baseline value
%                         index - peak index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Peak validation parameters
%             peakParam.  minPeakWidth - minimum peak width
%             ampThr - amplitude threshold; automatically determined if it
%             is zero

peakCount=length(peaks);
% Matrix pre-location
validatedPeaks=peaks;
minPeakWidth=peakParam.minPeakWidth;

if isempty(peakParam.ampThr)||peakParam.ampThr==false
    ampThr=getAmpThr(peaks,peakParam);
    peakParam.ampThr=ampThr;
else
    ampThr=peakParam.ampThr;
end

if peakParam.ampThr>0.9*max(SpSmooth)||peakParam.ampThr<1.1*min(SpSmooth)
    error('Peak validation threshold exceeds spectrum maximum and minimum values');
end

index=1;
for i=1:peakCount
    if peaks(i).endPos-peaks(i).startPos > minPeakWidth &&...
            peaks(i).maxVal-peaks(i).basl > ampThr
        validatedPeaks(index)=peaks(i);
        index=index+1;
    end
end

if index > 1;
    PeakCount=index-1;
    validatedPeaks=validatedPeaks(1:PeakCount);
else
    error('wrong peak peaking parameters: No Validated peaks')
end


minsegwidth=10.^10;;
for i=1:PeakCount
    startPos=validatedPeaks(i).startPos;
    maxPos=validatedPeaks(i).maxPos;
    endPos=validatedPeaks(i).endPos;
    segwidth=endPos-startPos;
    % Determine the peak boundaries
    edgeVal=validatedPeaks(i).maxVal.*peakParam.peakEdgeMax;
    LeftEdge=find(SpSmooth(startPos:maxPos)-validatedPeaks(i).basl>=...
        edgeVal,1,'first');
    if isempty(LeftEdge)
        validatedPeaks(i).LeftEdge=startPos;
    else
        validatedPeaks(i).LeftEdge=startPos+LeftEdge-1;
    end
    RightEdge=find(SpSmooth(maxPos:endPos)-validatedPeaks(i).basl>=...
        edgeVal,1,'last');
    if isempty(RightEdge)
        validatedPeaks(i).RightEdge=endPos;
    else
        validatedPeaks(i).RightEdge=maxPos+RightEdge-1;
    end
    if minsegwidth>segwidth
       minsegwidth=segwidth;
    end
end


if ~isempty(debug);
    peakStartPositions=[];
    peakEndPositions=[];
    peakMaxPositions=[];
    for i=1:PeakCount
        peakStartPositions=[peakStartPositions  validatedPeaks(i).LeftEdge];
        peakEndPositions=[peakEndPositions validatedPeaks(i).RightEdge];
        peakMaxPositions=[peakMaxPositions validatedPeaks(i).maxPos];
    end
%     plotPeaksOrSegments(SpSmooth,peakMaxPositions,...
%         peakStartPositions,peakEndPositions,debug)
end

return;
end

function pt = ppmToPt(ppmValues, firstPtPpm, resolution)

if(nargin < 2 || isempty(firstPtPpm))
    error('nargin < 2 || isempty(firstPtPpm)');
end;
if(nargin < 3 || isempty(resolution))
    resolution = ppmValues(2) - ppmValues(1);
end;

if(~isscalar(firstPtPpm))
    error('First ppm should be a number, got non-scalar value: %d', firstPtPpm);
end;
if(~isscalar(resolution))
    error('Resolution ppm should be a number, got non-scalar value: %d', resolution);
end;

ppmShift = ppmValues - firstPtPpm;

pt = round(ppmShift ./ resolution) + 1;

return
end


function [B,G] = sgolay(k,F,varargin)
%SGOLAY Savitzky-Golay Filter Design.
%   B = SGOLAY(K,F) designs a Savitzky-Golay (polynomial) FIR smoothing
%   filter B.  The polynomial order, K, must be less than the frame size,
%   F, and F must be odd.  
%
%   Note that if the polynomial order K equals F-1, no smoothing
%   will occur.
%
%   SGOLAY(K,F,W) specifies a weighting vector W with length F
%   containing real, positive valued weights employed during the
%   least-squares minimization.
%
%   [B,G] = SGOLAY(...) returns the matrix G of differentiation filters.
%   Each column of G is a differentiation filter for derivatives of order
%   P-1 where P is the column index.  Given a length F signal X, an
%   estimate of the P-th order derivative of its middle value can be found
%   from:
%
%                     ^(P)
%                     X((F+1)/2) = P!*G(:,P+1)'*X
%
%   See also SGOLAYFILT, FIR1, FIRLS, FILTER

%   References:
%     [1] Sophocles J. Orfanidis, INTRODUCTION TO SIGNAL PROCESSING,
%              Prentice-Hall, 1995, Chapter 8

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.12 $  $Date: 2002/04/15 01:17:19 $

error(nargchk(2,3,nargin));

% Check if the input arguments are valid
if round(F) ~= F, error('Frame length must be an integer.'), end
if rem(F,2) ~= 1, error('Frame length must be odd.'), end
if round(k) ~= k, error('Polynomial degree must be an integer.'), end
if k > F-1, error('The degree must be less than the frame length.'), end
if nargin < 3,
   % No weighting matrix, make W an identity
   W = eye(F);
else
   W = varargin{1};
   % Check for right length of W
   if length(W) ~= F, error('The weight vector must be of the same length as the frame length.'),end
   % Check to see if all elements are positive
   if min(W) <= 0, error('All the elements of the weight vector must be greater than zero.'), end
   % Diagonalize the vector to form the weighting matrix
   W = diag(W);
end

% Compute the projection matrix B
s = fliplr(vander(-(F-1)./2:(F-1)./2));
S = s(:,1:k+1);   % Compute the Vandermonde matrix

[Q,R] = qr(sqrt(W)*S,0);

G = S*inv(R)*inv(R)'; % Find the matrix of differentiators

B = G*S'*W; 


% [EOF] - sgolay.m
return;
end

function plotPeaksOrSegments(SpSmooth,peakMaxPositions,...
    peakStartPositions,peakEndPositions,debug,Marker)
% debugging purpose only

if nargin<6
    Marker=0.5;
end
% segmentCount=length(peakStartPositions);
% title(sprintf('%.d Segments Found ',segmentCount));
% hold on
% plot(SpSmooth,'b'); hold on;
% ylim = get(gca,'YLim');
% if ~isempty(peakMaxPositions)
%     line([peakMaxPositions;peakMaxPositions],...
%         [0.*ones([1 length(peakMaxPositions)]); ...
%         SpSmooth(peakMaxPositions)],'Linewidth', Marker,'Color','b');
% end
% line([peakStartPositions;peakStartPositions],...
%     [0.*ones([1 length(peakStartPositions)]);...
%     SpSmooth(peakStartPositions)*5],'Linewidth', Marker,'Color','g');
% line([peakEndPositions;peakEndPositions],...
%     [0.*ones([1 length(peakEndPositions)]);...
%     SpSmooth(peakEndPositions)*5],'Linewidth', Marker,'Color','k');
pause(debug);
hold off;
return
end


function x=simplepick(Sp,ppm,peakthresh)
next=1;
noise=std(Sp(1:500));
for i = 2:1:length(Sp)-1 %%find all peaks higher than given multiplier of noise%%
    if (Sp(i) > peakthresh*noise) && (Sp(i-1) < Sp(i)) && (Sp(i) >Sp(i+1))
        x(next)=i;
        next = next+1;
    end
end
end





