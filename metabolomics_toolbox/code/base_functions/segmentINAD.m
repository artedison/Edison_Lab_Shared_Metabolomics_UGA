function label=segmentINAD(X,XNoise,ppm1,ppm2)
% label=segment2D(X,XNoise,ppm1,ppm2,STNthresh,ppmbetweenNuc1,ppmbetweenNuc2)
%
% Spectral segmentation algorithm using Koradi noise calculation and
% morphological operations to grow segments - as outlined in HATS and
% mvaDANS
%
% Arguments:
% 
% X                  Data matrix of spectra
% XNoise             Data matrix of spectral noise from Koradi algorithm
% ppm1               Chemical shift vector of F2
% ppm2               Chemical shift vector of F1
% STNthresh          Signal to noise threshold for detection - default 10 for
%                    homonuclear, 4 for heteronuclear
% ppmbetweenNuc1     allowed distance between maxima in same crosspeak in F2 dimension, default 0.035 for homonuclear, 0.025 for heteronuclear
% ppmbetweenNuc2     allowed distance between maxima in same crosspeak in F1 dimension, default 0.035 for homonuclear, 0.3 for heteronuclear
%
% Return Values:
% label              matrix equal in size to spectrum with index values for each
%                    region 

% if abs(max(ppm1)-max(ppm2))<10
%     if exist('STNthresh')==0
%         STNthresh=10;
%     end
%     if exist('ppmbetweenNuc1')==0
%         ppmbetweenNuc1=0.035;
%     end
%     if exist('ppmbetweenNuc2')==0
%         ppmbetweenNuc2=0.035;
%     end
% else
%     if exist('STNthresh')==0
%         STNthresh=4;
%     end
%     if exist('ppmbetweenNuc1')==0
%         ppmbetweenNuc1=0.025;
%     end
%     if exist('ppmbetweenNuc2')==0
%         ppmbetweenNuc2=0.3;
%     end
% end
%

spect=mean(X,3);
noise=mean(XNoise,3);

%% changed by GSS for inadequate. commented out above code
STNthresh=20;
ppmbetweenNuc1=.3; %.05
ppmbetweenNuc2=.3;  %.5
%shift1=round(ppmbetweenNuc1/(ppm1(1)-ppm1(2)));
%shift2=round(ppmbetweenNuc2/(ppm2(1)-ppm2(2)));
shift1=4;
shift2=4;
%%
STN=spect./noise;

[FX,FY]=gradient(STN);
STNgradXY=sqrt(abs(FX.*FY));

prepeaksign1=zeros(size(spect));
for x=4:size(spect,1)-3
    for y=4:size(spect,2)-3
        if abs(ppm2(x)-ppm1(y))<=ppmbetweenNuc1
            continue
        elseif (spect(x,y)>=(max([spect(x-1,y-1),spect(x,y-1),spect(x+1,y-1),spect(x-1,y),spect(x+1,y),spect(x-1,y+1),spect(x,y+1),spect(x+1,y+1)])) && spect(x,y)>=abs(STNthresh*noise(x,y)));
            prepeaksign1(x,y)=1;
        end
    end
end


if isempty(ver('images'))==0
       
    se1=strel('rectangle',[2 2]);
    se2=strel('rectangle',[shift1 shift2]);
    
    peaksignbw=imdilate(prepeaksign1,strel(se2));   
    peaksign=bwlabel(peaksignbw); 
    
    %peaksign=corrchop(peaksign,prepeaksign1,X,shift1,shift2);
    
    %COMMENTED BY GSS FOR INADEQUATE
    %peaksign=chop_bins(peaksign,ppm1,ppm2);
    
    peaksign(peaksign==0)=Inf;
    
    for k=1:(shift1/2)
        peaksign=imerode(peaksign,se1);
        peaksign(isinf(peaksign)==1)=-Inf;
        peaksign=imdilate(peaksign,se1);
        peaksign(isinf(peaksign)==1)=Inf;
    end
    
    label=peaksign;
    label(isinf(label)==1)=0;
else
    factors1=divisor(size(X,1));
    factors2=divisor(size(X,2));
    [h1,k1]=min(abs(factors1-shift1));
    [h2,k2]=min(abs(factors2-shift2));
    shift1=factors1(k1);
    shift2=factors2(k2);
    
    peaksignbw=maxfilt2(prepeaksign1,[shift1,shift2],'same');   
    peaksign=slrlabel(peaksignbw);
    
    %peaksign=corrchop(peaksign,prepeaksign1,X,shift1,shift2);
    peaksign=chop_bins(peaksign,ppm1,ppm2);
    
    peaksign(peaksign==0)=Inf;
    
    for k=1:(shift1/2)
        peaksign=minfilt2(peaksign);
        peaksign(isinf(peaksign)==1)=-Inf;
        peaksign=maxfilt2(peaksign);
        peaksign(isinf(peaksign)==1)=Inf;
    end
    
    label=peaksign;
    label(isinf(label)==1)=0;
end
return 
thresh=3;
levels=10;
range=3;
[FXlab,FYlab]=gradient(label);
figure, imagesc(ppm1,ppm2,FXlab+FYlab)
cmapregions=[0.847058832645416,0.160784319043159,0;1,1,1;0.847058832645416,0.160784319043159,0];
colormap(cmapregions);
vector=(2.^[-1*range:(range-(-1*range))/(levels-1):range])*(thresh*std(std(spect./noise)));
caxis([-1 1])
hold on
h=contour(ppm1,ppm2,spect./noise,vector,'EdgeColor','k');
set(gca,'XDir','rev')
set(gca,'YDir','rev')
end


function label=chop_bins(label,ppm1,ppm2,dirmax)

maxlabel=max(max(label))+1;
loopcount=1;

if exist('dirmax')==0
if abs(max(ppm1)-max(ppm2))<10
    indirmax=.15;
    dirmax=.15;
else
    indirmax=10; %1
    dirmax=20; %2 %changed by gss for inadequate from 1
end
else
    indirmax=dirmax;
end

while loopcount<=max(max(label))
    
    [I,J]=find(label==loopcount);
     
    maxdir=max(ppm1(J));
    mindir=min(ppm1(J));
    maxindir=max(ppm2(I));
    minindir=min(ppm2(I));
    dirsize=maxdir-mindir;
    indirsize=maxindir-minindir;
    
    if dirsize>=dirmax && dirsize>indirsize
        relabel=find(ppm1(J)>(maxdir-(dirsize)/2));
        relabelsubs=sub2ind(size(label),I(relabel),J(relabel));
        label(relabelsubs)=maxlabel;
        maxlabel=maxlabel+1;
        disp(loopcount)
    elseif indirsize>=indirmax && indirsize>dirsize
        relabel=find(ppm2(I)>(maxindir-(indirsize)/2));
        relabelsubs=sub2ind(size(label),I(relabel),J(relabel));
        label(relabelsubs)=maxlabel;
        maxlabel=maxlabel+1;
    end
    loopcount=loopcount+1;
end

end


function d = divisor(n)
%% divisor : provides a list of integer divisors of a number.
% divisor(n) : row vector of all distinct divisors of a positive integer N, 
%               including 1 and N.
%
% Remark:
%   This function uses the default factor() routine in Matlab and hence is 
% limited to input values upto 2^32. However if factor() routine does get
% updated for larger integers, this function will still work fine.
%   Using factor() provides a significant speed improvement over manually 
% seaching for the each divisor of n.
%
% Example:
%   a = divisor(12);
%   returns -> a = [1, 2, 3, 4, 6, 12];
%
% See Also:
%   factor, primes

% Author: Yash Kochar ( yashkochar@yahoo.com )
% Last modified: 21st June 2009
%-------------------------------------------------------------------------------

%% Input error check :
%   Check whether input is positive integer and scalar.
if ~isscalar(n)
    error('divisor:NonScalarInput','Input must be a scalar.');
end
if (n < 1) || (floor(n) ~= n)
  error('divisor:PositiveIntegerOnly', 'Input must be a positive integer.'); 
end

%% Find prime factors of number :
pf = factor(n);         % Prime factors of n
upf = unique(pf);       % Unique

%% Calculate the divisors
d = upf(1).^(0:1:sum(pf == upf(1)))';
for f = upf(2:end)
    d = d*(f.^(0:1:sum(pf == f)));
    d = d(:);
end
d = sort(d)';   % To further improve the speed one may remove this sort command
                %   Just remember to take the transpose of "d" to get a result
                %   as a row vector instead of a column vector.
end



function Y = maxfilt2(X,varargin)
%  MAXFILT2    Two-dimensional max filter
%
%     Y = MAXFILT2(X,[M N]) performs two-dimensional maximum
%     filtering on the image X using an M-by-N window. The result
%     Y contains the maximun value in the M-by-N neighborhood around
%     each pixel in the original image. 
%     This function uses the van Herk algorithm for max filters.
%
%     Y = MAXFILT2(X,M) is the same as Y = MAXFILT2(X,[M M])
%
%     Y = MAXFILT2(X) uses a 3-by-3 neighborhood.
%
%     Y = MAXFILT2(..., 'shape') returns a subsection of the 2D
%     filtering specified by 'shape' :
%        'full'  - Returns the full filtering result,
%        'same'  - (default) Returns the central filter area that is the
%                   same size as X,
%        'valid' - Returns only the area where no filter elements are outside
%                  the image.
%
%     See also : MINFILT2, VANHERK
%

% Initialization
[S, shape] = parse_inputs(varargin{:});

% filtering
Y = vanherk(X,S(1),'max',shape);
Y = vanherk(Y,S(2),'max','col',shape);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S, shape] = parse_inputs(varargin)
shape = 'same';
flag = [0 0]; % size shape

for i = 1 : nargin
   t = varargin{i};
   if strcmp(t,'full') & flag(2) == 0
      shape = 'full';
      flag(2) = 1;
   elseif strcmp(t,'same') & flag(2) == 0
      shape = 'same';
      flag(2) = 1;
   elseif strcmp(t,'valid') & flag(2) == 0
      shape = 'valid';
      flag(2) = 1;
   elseif flag(1) == 0
      S = t;
      flag(1) = 1;
   else
      error(['Too many / Unkown parameter : ' t ])
   end
end

if flag(1) == 0
   S = [3 3];
end
if length(S) == 1;
   S(2) = S(1);
end
if length(S) ~= 2
   error('Wrong window size parameter.')
end

end
end



function Y = minfilt2(X,varargin)
%  MINFILT2    Two-dimensional min filter
%
%     Y = MINFILT2(X,[M N]) performs two-dimensional minimum
%     filtering on the image X using an M-by-N window. The result
%     Y contains the minimun value in the M-by-N neighborhood around
%     each pixel in the original image. 
%     This function uses the van Herk algorithm for min filters.
%
%     Y = MINFILT2(X,M) is the same as Y = MINFILT2(X,[M M])
%
%     Y = MINFILT2(X) uses a 3-by-3 neighborhood.
%
%     Y = MINFILT2(..., 'shape') returns a subsection of the 2D
%     filtering specified by 'shape' :
%        'full'  - Returns the full filtering result,
%        'same'  - (default) Returns the central filter area that is the
%                   same size as X,
%        'valid' - Returns only the area where no filter elements are outside
%                  the image.
%
%     See also : MAXFILT2, VANHERK
%

% Initialization
[S, shape] = parse_inputs(varargin{:});

% filtering
Y = vanherk(X,S(1),'min',shape);
Y = vanherk(Y,S(2),'min','col',shape);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S, shape] = parse_inputs(varargin)
shape = 'same';
flag = [0 0]; % size shape

for i = 1 : nargin
   t = varargin{i};
   if strcmp(t,'full') & flag(2) == 0
      shape = 'full';
      flag(2) = 1;
   elseif strcmp(t,'same') & flag(2) == 0
      shape = 'same';
      flag(2) = 1;
   elseif strcmp(t,'valid') & flag(2) == 0
      shape = 'valid';
      flag(2) = 1;
   elseif flag(1) == 0
      S = t;
      flag(1) = 1;
   else
      error(['Too many / Unkown parameter : ' t ])
   end
end

if flag(1) == 0
   S = [3 3];
end
if length(S) == 1;
   S(2) = S(1);
end
if length(S) ~= 2
   error('Wrong window size parameter.')
end

end
end

function Y = vanherk(X,N,TYPE,varargin)
%  VANHERK    Fast max/min 1D filter
%
%    Y = VANHERK(X,N,TYPE) performs the 1D max/min filtering of the row
%    vector X using a N-length filter.
%    The filtering type is defined by TYPE = 'max' or 'min'. This function
%    uses the van Herk algorithm for min/max filters that demands only 3
%    min/max calculations per element, independently of the filter size.
%
%    If X is a 2D matrix, each row will be filtered separately.
%    
%    Y = VANHERK(...,'col') performs the filtering on the columns of X.
%    
%    Y = VANHERK(...,'shape') returns the subset of the filtering specified
%    by 'shape' :
%        'full'  - Returns the full filtering result,
%        'same'  - (default) Returns the central filter area that is the
%                   same size as X,
%        'valid' - Returns only the area where no filter elements are outside
%                  the image.
%
%    X can be uint8 or double. If X is uint8 the processing is quite faster, so
%    dont't use X as double, unless it is really necessary.
%

% Initialization
[direc, shape] = parse_inputs(varargin{:});
if strcmp(direc,'col')
   X = X';
end
if strcmp(TYPE,'max')
   maxfilt = 1;
elseif strcmp(TYPE,'min')
   maxfilt = 0;
else
   error([ 'TYPE must be ' char(39) 'max' char(39) ' or ' char(39) 'min' char(39) '.'])
end

% Correcting X size
fixsize = 0;
addel = 0;
if mod(size(X,2),N) ~= 0
   fixsize = 1;
   addel = N-mod(size(X,2),N);
   if maxfilt
      f = [ X zeros(size(X,1), addel) ];
   else
      f = [X repmat(X(:,end),1,addel)];
   end
else
   f = X;
end
lf = size(f,2);
lx = size(X,2);
clear X

% Declaring aux. mat.
g = f;
h = g;

% Filling g & h (aux. mat.)
ig = 1:N:size(f,2);
ih = ig + N - 1;

g(:,ig) = f(:,ig);
h(:,ih) = f(:,ih);

if maxfilt
   for i = 2 : N
      igold = ig;
      ihold = ih;
      
      ig = ig + 1;
      ih = ih - 1;
      
      g(:,ig) = max(f(:,ig),g(:,igold));
      h(:,ih) = max(f(:,ih),h(:,ihold));
   end
else
   for i = 2 : N
      igold = ig;
      ihold = ih;
      
      ig = ig + 1;
      ih = ih - 1;
      
      g(:,ig) = min(f(:,ig),g(:,igold));
      h(:,ih) = min(f(:,ih),h(:,ihold));
   end
end
clear f

% Comparing g & h
if strcmp(shape,'full')
   ig = [ N : 1 : lf ];
   ih = [ 1 : 1 : lf-N+1 ];
   if fixsize
      if maxfilt
         Y = [ g(:,1:N-1)  max(g(:,ig), h(:,ih))  h(:,end-N+2:end-addel) ];
      else
         Y = [ g(:,1:N-1)  min(g(:,ig), h(:,ih))  h(:,end-N+2:end-addel) ];
      end
   else
      if maxfilt
         Y = [ g(:,1:N-1)  max(g(:,ig), h(:,ih))  h(:,end-N+2:end) ];
      else
         Y = [ g(:,1:N-1)  min(g(:,ig), h(:,ih))  h(:,end-N+2:end) ];
      end
   end
   
elseif strcmp(shape,'same')
   if fixsize
      if addel > (N-1)/2
%          disp('hoi')
         ig = [ N : 1 : lf - addel + floor((N-1)/2) ];
         ih = [ 1 : 1 : lf-N+1 - addel + floor((N-1)/2)];
         if maxfilt
            Y = [ g(:,1+ceil((N-1)/2):N-1)  max(g(:,ig), h(:,ih)) ];
         else
            Y = [ g(:,1+ceil((N-1)/2):N-1)  min(g(:,ig), h(:,ih)) ];
         end
      else   
         ig = [ N : 1 : lf ];
         ih = [ 1 : 1 : lf-N+1 ];
         if maxfilt
            Y = [ g(:,1+ceil((N-1)/2):N-1)  max(g(:,ig), h(:,ih))  h(:,lf-N+2:lf-N+1+floor((N-1)/2)-addel) ];
         else
            Y = [ g(:,1+ceil((N-1)/2):N-1)  min(g(:,ig), h(:,ih))  h(:,lf-N+2:lf-N+1+floor((N-1)/2)-addel) ];
         end
      end            
   else % not fixsize (addel=0, lf=lx) 
      ig = [ N : 1 : lx ];
      ih = [ 1 : 1 : lx-N+1 ];
      if maxfilt
         Y = [  g(:,N-ceil((N-1)/2):N-1) max( g(:,ig), h(:,ih) )  h(:,lx-N+2:lx-N+1+floor((N-1)/2)) ];
      else
         Y = [  g(:,N-ceil((N-1)/2):N-1) min( g(:,ig), h(:,ih) )  h(:,lx-N+2:lx-N+1+floor((N-1)/2)) ];
      end
   end      
   
elseif strcmp(shape,'valid')
   ig = [ N : 1 : lx];
   ih = [ 1 : 1: lx-N+1];
   if maxfilt
      Y = [ max( g(:,ig), h(:,ih) ) ];
   else
      Y = [ min( g(:,ig), h(:,ih) ) ];
   end
end

if strcmp(direc,'col')
   Y = Y';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [direc, shape] = parse_inputs(varargin)
direc = 'lin';
shape = 'same';
flag = [0 0]; % [dir shape]

for i = 1 : nargin
   t = varargin{i};
   if strcmp(t,'col') & flag(1) == 0
      direc = 'col';
      flag(1) = 1;
   elseif strcmp(t,'full') & flag(2) == 0
      shape = 'full';
      flag(2) = 1;
   elseif strcmp(t,'same') & flag(2) == 0
      shape = 'same';
      flag(2) = 1;
   elseif strcmp(t,'valid') & flag(2) == 0
      shape = 'valid';
      flag(2) = 1;
   else
      error(['Too many / Unkown parameter : ' t ])
   end
end
end
end


function label=slrlabel(im,connectivity)

if exist('connectivity')==0
    connectivity=4;
end

[h(:,1),h(:,2)]=find(im==1);
index=sub2ind(size(im),h(:,1),h(:,2));

% tic
if connectivity==4
pairs=zeros(4*size(h,1),2);
else
    pairs=zeros(8*size(h,1),2);
end
z=1;
for k=1:size(h,1)
    if connectivity==8
    connected=intersect(find(abs(h(k,1)-h(:,1))<=1),find(abs(h(k,2)-h(:,2))<=1));
    elseif connectivity==4
    connected=union(intersect(find(abs(h(k,1)-h(:,1))<=0),find(abs(h(k,2)-h(:,2))<=0)),union(intersect(find(abs(h(k,1)-h(:,1))<=0),find(abs(h(k,2)-h(:,2))<=1)),intersect(find(abs(h(k,1)-h(:,1))<=1),find(abs(h(k,2)-h(:,2))<=0))));
    else
            connected=union(intersect(find(abs(h(k,1)-h(:,1))<=0),find(abs(h(k,2)-h(:,2))<=1)),intersect(find(abs(h(k,1)-h(:,1))<=1),find(abs(h(k,2)-h(:,2))<=0)));

    end
    
    for w=1:length(connected)
        pairs(z,:)=[k,connected(w)];
        z=z+1;
    end
    clear connected
end
pairs(find(pairs(:,1)==0),:)=[];

% toc

if size(pairs,1)<4000
    A = accumarray(pairs, 1);
else
    A = sparse(pairs(:,1), pairs(:,2), 1);
end

[p,q,r,s] = dmperm(A);

label=zeros(size(im));

for k=1:(size(r,2)-1)
    indices=index(p(r(k):r(k+1)-1));
    label(indices)=k;
    clear indices
end
end

