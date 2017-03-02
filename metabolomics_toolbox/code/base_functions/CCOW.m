function XWarped = CCOW(X, T, varargin)
% Constrained Correlation Optimized Warping for aligning 1D-NMR spectra of
% complex biological\chemical mixtures

% Author: Kirill Veselkov, Imperial College 2009; The code is adapted and modified for efficiency from
% Giorgio Tomasi / Frans van den Berg (Quality and Technology, Department of Food Science, Faculty of Life Sciences, University of Copenhagen) 

% References: Niels-Peter Vest Nielsen, Jens Micheal Carstensen and Jørn Smedegaard 'Aligning of singel and multiple
%             wavelength chromatographic profiles for chemometric data analysis using correlation optimised warping' 
%             J. Chrom. A 805(1998)17-35
%
%             Correlation optimized warping and dynamic time warping as preprocessing methods for chromatographic Data
%             Giorgio Tomasi, Frans van den Berg and Claus Andersson, Journal of Chemometrics 18(2004)231-241
%

% Input: X [nSmpls x nVrbls] - 1D-NMR spectra of biological samples
%        T [1 x nVrbls] - the target spectrum with respect to peak positions of
%        which peak positions of other samples are to be corrected
%        'maxPeakShift' - max peak position adjustment
%        'SegWidth'     - the segment width
% Output: XWarped - warped spectra of biological samples

[maxPeakShift,SegLength,CCpower,Slack]= setPrmtrsCCOW(varargin);

if nargin<2
    error('incorrect number of input arguments');
end

%% Initialize the algorithm parameters:
[nSmpls,nVrbls]      = size(X);                  % nSmpls - the number of samples
                                                 % nVrbls - the number of variables
nSgnts               = floor((nVrbls-1)/(SegLength-1));  % nSgnts - the number of segments into which spectra are divided
if nSgnts < 2
    error('ERROR: The segment length is larger than the number of spectral data points');
end
SegLengths(1:nSgnts) = SegLength-1;
SegLengths(nSgnts)   = SegLengths(nSgnts)+rem(nVrbls-1,SegLength-1); %add remaining points to the last segment
SegBndrs = cumsum([1,SegLengths]);                               %initial segment boundaries

if maxPeakShift==0
    error('ERROR: The maximum segment shift must not be set to zero');
    return;
end
XWarped         = repmat(NaN,[nSmpls,nVrbls]);
WarpedBndrPstns = repmat(NaN,[nSmpls,nSgnts+1]);

%% Calculate the warping constraints on segment boundaries:
if isempty(Slack)
    if any(maxPeakShift + 2 > SegLengths)
        Slack = SegLength - 4; % the slack parameter is set to the segment length -4
    else
        Slack = maxPeakShift;
    end
elseif any(Slack +2 > SegLengths)
    Slack = SegLength - 3;
    warning('The slack parameter must be smaller that the average peak width by at least two data points');
end

offs              = (Slack * [-1,1]') * (0:nSgnts);                             % Calculate the COW constraints on segment boundaries. Under these constraints,
WarpCntrsStart    = SegBndrs(ones(2,1),1:nSgnts+1) + offs;                      % the flexibility of segment boundaries (and thus the peak position correction)
WarpCntrsEnd      = SegBndrs(ones(2,1),1:nSgnts+1) + offs(:,nSgnts + 1:-1:1);   % increases towards the middle of the chemical shift scale. This is not of a
WarpCntrsCOW(1,:) = max(WarpCntrsStart(1,:),WarpCntrsEnd(1,:));                 % a particular advantage for NMR biological spectra in which the extent of variable 
WarpCntrsCOW(2,:) = min(WarpCntrsStart(2,:),WarpCntrsEnd(2,:));                 % peak positions is roughly similar across the chemical shift scale.

if isempty(maxPeakShift)
    WarpCntrs(1,:) = WarpCntrsCOW(1,:);
    WarpCntrs(2,:) = WarpCntrsCOW(2,:);
else
    WarpCntrs(1,:) = max(WarpCntrsCOW(1,:),SegBndrs - maxPeakShift); % Apply additional constraints on segment boundaries to maintain comparable adjustment of peak positions
    WarpCntrs(2,:) = min(WarpCntrsCOW(2,:),SegBndrs + maxPeakShift); % across the chemical shift scale
end

%% Initialize a table for the dynamic programming:
nBndrPsnts = cumsum([0,diff(WarpCntrs) + 1]);
DPTable = repmat(0,[nBndrPsnts(end),nSmpls,2]); % Table: each column relates to the boundary of a segment
                                                %  (:,i,1) the value of the cumulative benefit function of the boundary
                                                %  of the i-th segment
                                                %  (:,i,2) pointer to the index of the preceding optimal boundary point

DPTable(2:end,:,1) = -Inf;                      % All values of the cumulative benefit function are initially set to -Inf
SegBndrPstns = repmat(1,[1,nBndrPsnts(end)]);   % Pre-allocate a vector for the storage of segment boundary positions

%% Calculate parameters for the piece-wise linear interpolation:
DiffXXi = cell(nSgnts,1);
XIndices = DiffXXi;
Xdiff = diff(X,1,2);
Xdiff(:,size(Xdiff,2)+1) = 0;
uniqSegLengths  = unique(SegLengths);
segSlacks         = -Slack:Slack;
for uniqSegLength = uniqSegLengths
    segIndcs                 = find(SegLengths==uniqSegLength);
    [iSegDiffXXi,XIndicesSeg] = Interpolate(uniqSegLength+1,segSlacks);
    [DiffXXi{segIndcs}]      = deal(iSegDiffXXi);
    [XIndices{segIndcs}]     = deal(XIndicesSeg);
end
warning('off','MATLAB:divideByZero')              % Turn off the division by zero warning

%% Initialize the forward phase of the dynamic programming algorithm:
for iSeg=1:nSgnts
    iSegWarps                 = segSlacks + SegLengths(iSeg);            % all possible warpings of the i-th segment with the given length
    iSSegBndrPstns            = WarpCntrs(1,iSeg+1):WarpCntrs(2,iSeg+1); % get indices of all boundary positions of the i+1 segment
    iStart                    = nBndrPsnts(iSeg+1)+1;
    iEnd                      = nBndrPsnts(iSeg+2);
    SegBndrPstns(iStart:iEnd) = iSSegBndrPstns;                       % store all segment boundary positions
    iSegXIndices = XIndices{iSeg}';                                   % Xi indices used in the interpolation of the i-th segment 
    iSegDiffXXi = DiffXXi{iSeg}';                                     % Difference between X and Xi used in the interpolation of the i-th segment
        
    tSegX      = T(SegBndrs(iSeg):SegBndrs(iSeg+1));                  % the intensity values of the target segment
    tCntrdSeg  = tSegX - mean(tSegX);                                 % the centered intensity values of the targe segment used in the calculation of the correlation coefficients
    tCntrdSegNorm  = norm(tCntrdSeg);                                 % the norm of the target segment used in the denominator of the correlation coefficient
    iBndPstnCount=iStart;
    for iSSegBndrPstn = iSSegBndrPstns
        iSegBndrPstns   = iSSegBndrPstn - iSegWarps;                 % the possible boundary positions of the i-th segment 
        iSegAccBndPstns = iSegBndrPstns >= WarpCntrs(1,iSeg) & iSegBndrPstns <= WarpCntrs(2,iSeg); % the accepted boundary positions of the i-th segment
        PrevBndrIndices = nBndrPsnts(iSeg) + 1 + iSegBndrPstns(iSegAccBndPstns)-WarpCntrs(1,iSeg);  % store index of the previous segment in a table
        niSegAccBndPstns    = sum(iSegAccBndPstns);                   % the number of acceptable boundary positions for the i-th segment
        if (niSegAccBndPstns)
         iSegInterXIndices             = iSSegBndrPstn + iSegXIndices(:,iSegAccBndPstns); % Interpolation signal indexes for all the allowed arcs for node i_node
         iSegDiffAccXXi                = iSegDiffXXi(:,iSegAccBndPstns);    % Interpolation coefficients for all the allowed arcs for node i_node
         iSegDiffAccXXi                = iSegDiffAccXXi(:)';
         iSegDiffAccXXi                = iSegDiffAccXXi(ones(nSmpls,1),:);
         iSegX                         = X(:,iSegInterXIndices);
         iSegXdiff                     = Xdiff(:,iSegInterXIndices);
         iSegX                         = reshape((iSegX + iSegDiffAccXXi .* iSegXdiff)',SegLengths(iSeg)+1,niSegAccBndPstns * nSmpls);       % Interpolate for all allowed predecessors
         iSegMean                      = sum(iSegX)/size(iSegX,1); 
         iCntrdSegNorm                 = sqrt(sum(iSegX.^2) - size(iSegX,1) * iSegMean.^2);      
         iSegsCCs                      = (tCntrdSeg * iSegX)./(tCntrdSegNorm * iCntrdSegNorm)+1;   
         iSegsCCs(~isfinite(iSegsCCs)) = 0;
         iSegsCCs                      = reshape(iSegsCCs,niSegAccBndPstns,nSmpls);
         CostFun                       = DPTable(PrevBndrIndices,:,1) + iSegsCCs.^CCpower;   % store the optimal value of loss function from all predecessors
         [CostFunValue,index]          = max(CostFun,[],1);                                  % optimal value of loss function from all predecessors
         DPTable(iBndPstnCount,:,1)    = CostFunValue;
         DPTable(iBndPstnCount,:,2)    = PrevBndrIndices(index);
         iBndPstnCount                 = iBndPstnCount+1;
        else
            error('consider the case');
        end
    end
end

%% Initialize the backward phase of the dynamic programming algorithm:
for iSmpl = 1:nSmpls                                  
   OptPredSegBndrIndex           = length(SegBndrPstns);
   WarpedBndrPstns(iSmpl,nSgnts + 1) = nVrbls;
   for iSegBndr = nSgnts:-1:1
      OptPredSegBndrIndex             = DPTable(OptPredSegBndrIndex,iSmpl,2);
      WarpedBndrPstns(iSmpl,iSegBndr) = SegBndrPstns(OptPredSegBndrIndex);
      ntSegXiPnts                     = SegBndrs(iSegBndr+1)-SegBndrs(iSegBndr)+1;
      iSegWarpIndcs                   = WarpedBndrPstns(iSmpl,iSegBndr):WarpedBndrPstns(iSmpl,iSegBndr+1);
      niSegWarpPnts                   = WarpedBndrPstns(iSmpl,iSegBndr+1) - WarpedBndrPstns(iSmpl,iSegBndr)+1;
      iSegXi                          = 1:niSegWarpPnts;
      iSegYi                          = X(iSmpl,iSegWarpIndcs);
      iSegX                           =(0:ntSegXiPnts-1) * (niSegWarpPnts-1)/(ntSegXiPnts-1) + 1;
      iSegY                           = interp1q(iSegXi',iSegYi',iSegX');
      XWarped(iSmpl,SegBndrs(iSegBndr):SegBndrs(iSegBndr+1))        = iSegY;
   end
end

return;

%% Calculate parameters for the piece-wise linear interpolation
function [DiffXXi,XIndices] = Interpolate(nSegPnts,segSlacks)
nSS          = length(segSlacks);
segWarpLngth = nSegPnts + segSlacks;
DiffXXi      = zeros(nSS,nSegPnts); % difference between X and Xi used in the interpolation
XIndices     = zeros(nSS,nSegPnts); 

for iSS = 1:nSS
X                          = (0:nSegPnts-1) * (segWarpLngth(iSS)-1)/(nSegPnts-1) + 1;
Xi                         = floor(X);
Xi(X >= segWarpLngth(iSS)) = segWarpLngth(iSS) - 1;
DiffXXi(iSS,:)             = X - Xi;
XIndices(iSS,:)            = Xi - segSlacks(iSS)-nSegPnts; % Xi-1 = Xi - segLength - Slack
end
return

function [maxPeakShift,SegLength,CCpower,Slack]= setPrmtrsCCOW(inputPrmtrs)
% default parameters
maxPeakShift = 40;
SegLength    = 25;
CCpower      = 1;
Slack        = []; 
nPrmtrs=length(inputPrmtrs);
for prmtrInd=1:2:nPrmtrs
    if  strcmp('maxPeakShift',inputPrmtrs{prmtrInd})
        maxPeakShift=inputPrmtrs{prmtrInd+1};      
    elseif strcmp('Slack',inputPrmtrs{prmtrInd})
        Slack = inputPrmtrs{prmtrInd+1};
    elseif strcmp('SegLength',inputPrmtrs{prmtrInd})
        SegLength=inputPrmtrs{prmtrInd+1};
    elseif  strcmp('CCpower',inputPrmtrs{prmtrInd})
        CCpower=inputPrmtrs{prmtrInd+1};      
    end
end

return;