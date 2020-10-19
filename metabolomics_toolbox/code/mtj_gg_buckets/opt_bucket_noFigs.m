function [ZNN,Z,I_b,S_b]=opt_bucket_noFigs(ppm,X,size_bucket,slackness)
%% opt_bucket_noFigs
%
    % Author: SAA Sousa (MTJ modified)
    % Version: 0.1
    % Tested on Matlab Version R2020b
    % Date: JUL2020
    %
    % Description:
    %       Version of opt_bucket, modified to skip figure generation for
    %       speed increase within optimize_optBucket(). Modification
    %       entailed commenting out the call to internal function 
    %       visualize() ("visualize(ZNN, I_b,ppm,X)" on Line 194. 
    %
    % Input:
    %       see original documentation below
    %
    % Output:
    %       see original documentation below
    %
    % Log:
    %
    % Example run:
    %
    %       [ZNN,Z,I_b,S_b]=opt_bucket_noFigs(ppm,X,size_bucket,slackness);
    %
%% Original Documentation: 
%
%
% size_bucket = .01
% slackness = .2
% ====================================================================================
%     Optmized Bucketing Algorithm (OBA) for bucketing in NMR 1H spectra defining local minima as buckets boundaries.
%     Usage: [ZNN,Z,I_b,S_b]=opt_bucket(xaxis,X,size_bucket,slackness);
% ====================================================================================
% ============================================================================
%   Copyright (c) 2013 Laboratory for Theoretical and Applied Chemometrics
%    Institute of Chemistry - University of Campinas
%	Please visit http://lqta.iqm.unicamp.br/
%  	If you publish data that made use of this method, you are also required to
%  	show this copyright notice. You are free to extend or modify this algorithm, however you are required
%  	to keep this copyright notice in your projects. Commercial use is explicitly prohibited, contact us for commercial options
%      ============================================================================
% ====================================================================================
% inputs:
%       ppm         = vector with chemical shifts in ppm;
%       X           = Matrix with row-vectors relatives to each NMR spectrum of each sample;
%       size_bucket = initial size to each bucket (e.g., 0.05 ppm);
%       slackness   = percentage of the size_bucket, it ranges from 0 to 1 (e.g.,slackness=0.5 or 50%). 
%                     If you set slackness=0, the algorithm performs the "conventional bucketing".
%     
% outputs:
%       ZNN         = Matrix with non-normalized buckets. rows = samples columns = buckets;
%       Z           = Matrix with normalized integrals to one (areas) of each bucket. rows
%                   = samples, columns = buckets;
%       I_b         = Chemical shift intervals to each bucket;
%       S_b         = Final size of each bucket;
%
% @Author: SAA Sousa // LQTA - UNICAMP // (c) 2013 // www.lqta.iqm.unicamp.br
% To contact the author, please e-mail to samuelquimica@gmail.com
% @Reference:  	SAA Sousa, A Magalhães, MMC Ferreira, Optimized Bucketing for NMR spectra: Three
%               case studies. Chemometrics and Intelligent Laboratory Systems, 122, 93-102 (2013). 

%% ====================================================================================
% Copyright Notice
% This is a license to use and modify SOFTWARE 
% produced by Laboratory for Theoretical and Applied Chemometrics
%   - Institute of Chemistry - University of Campinas
% If you use, modify or redistribute the SOFTWARE and/or
% your source modifications, you agree:
% (1) to use the SOFTWARE and/or your source modifications
% solely as part of your research and not in any commercial
% product;
% (2) that the SOFTWARE and/or your source modifications
% will not be distributed for profit;
% (3) all copyright notices and this license note with the
% SOFTWARE are retained any redistribution of
% the SOFTWARE or any portion thereof;
% (4) to indemnify, hold harmless, and defend the Laboratory for Theoretical 
%     and Applied Chemometrics- Institute of Chemistry - University of Campinas
% from and against any claims or lawsuits that arise or
% result from the use of the SOFTWARE or your source
% modifications;
% (5) to properly reference the SOFTWARE when used in
% your research in any publication that may result from
% that research.
%
% Reserved Rights. The Laboratory for Theoretical and Applied Chemometrics
% - Institute of Chemistry - University of Campinas retain title and all ownership
% rights to the SOFTWARE & DATA.
%
% Copyright 2013 Laboratory for Theoretical and Applied Chemometrics
%   - Institute of Chemistry - University of Campinas
% @Reference:  	SAA Sousa, A Magalhães, MMC Ferreira, Optimized Bucketing for NMR spectra: Three
%               case studies. Chemometrics and Intelligent Laboratory
%               Systems, 122, 93-102 (2013). 
%%
% ====================================================================================
% Step 0: Checking the inputs. The program only works with four inputs.
% =========================================================================
if (nargin < 4)
    help opt_bucket;
    return;
end
reverse = 0;
if ppm(1)-ppm(2)<0
    reverse = 1;
    ppm=fliplr(ppm);
    X=fliplr(X);
end


% ====================================================================================
% Step 1: Reading the X dimensions (Matrix with NMR spectra);
% =========================================================================
[p,q]=size(X);
% ====================================================================================
% Step 2: b = interval sampling;
% =========================================================================
b = ppm(1)-ppm(2);
% ====================================================================================
% Step 3: a = number of intervals to the buckets;
% =========================================================================
a = size_bucket./b;
a = round(a);
% ====================================================================================
% Step 4: l = slackness in number of intervals;
% =========================================================================
l=slackness*a;
l=round(l);
% ====================================================================================
% Step 5: R = Averaging NMR spectrum of the all samples;
% =========================================================================
R=mean(X);
% ====================================================================================
% Step 6: local minima definition;
% =========================================================================
v=[];
for t=1+a:a:q-a;
    [~,I]=min(R(t-l:t+l));
    f=((t-(1+a))./a)+1;
    v(f)=I+a*(f-1)+(a-l);
end
% ====================================================================================
% Step 7: removing intersections and calculating the integrals (areas);
% =========================================================================
z = unique(v);
v = [1 z q];
Z = zeros(p,(length(v)-1));
for j = 1:p;
    for n = 1:(length(v)-1);
        x_prov = X(j,v(n):v(n+1));
        Z(j,n) = trapzeq(x_prov);
    end
end
vv = zeros(p,(length(v)));
for j = 1:p;
    vv(j,:) = X(j,v(:));
end
for tt = 2:(length(v)-1);
    Z(:,tt) = Z(:,tt) - vv(:,tt);
end
ZNN=Z;
% ====================================================================================
% Step 8: Normalize to unit area;
% =========================================================================
m=size(Z,1);
for k=1:m;
    Z(k,:)=Z(k,:)./(sum(Z(k,:)));
end
% ====================================================================================
% Step 9: Reading the optimized boundaries;
% =========================================================================
A=[];
for k=1:length(v);
    A(k)=ppm(v(k));
end
s=length(A);
I_b=[A(1:(s-1))',A(2:s)'];
% ====================================================================================
% Step 10: Final bucket size;
% =========================================================================
T=[];
for k=1:(s-1);
    T(k)=A(k)-A(k+1);
    S_b=T';
end
if reverse
    ZNN=fliplr(ZNN);
    Z=fliplr(Z);
    I_b=fliplr(I_b);
    S_b=fliplr(S_b);
end
%visualize(ZNN, I_b,ppm,X)
end
%% ************************************************************************
% ====================================================================================
% helper function;
% =========================================================================
function q = trapzeq(y)
n = length(y);
sums = y(1) + 2*sum(y(2:n-1)) + y(n);
q = sums/2;
end
%% ************************************************************************

function visualize(ZNN, I_b,ppm,X)
plotr(fliplr(mean(I_b,2)'),ZNN)
h1=gca;
figure, plotr(ppm,X)
h2=gca;
for i=1:length(I_b)
    line([I_b(i,1), I_b(i,1)],[0 max(max(X))])
end
linkaxes([h1,h2],'x')
end
