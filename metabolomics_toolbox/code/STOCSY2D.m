function S = STOCSY2D(H, X, ppmH, ppmX, loadingsH, loadingsX, varargin)
% *********************
% STOCSY2D
% Insert function description here
%
% Required arguments:
%     H:
%     X:
%     ppmH:
%     ppmX:
% Optional arguments:
%     loadingsH:
%     loadingsX:
% Optional Name/Value Pairs:
%     overlayHSQC: path to HSQC .ft file
%     threshH:
%     threshX:
%     numPointsH:
%     numPointsX:
%     threshCorr:
%     corrMode: 'pos', 'neg', 'both' (default)
%     peakPickX: 1 or 0 (default)
%     peakPickH: 1 or 0 (default)
%
% *********************
% Examples:
% STOCSY2D(H, X, ppmH, ppmX, 'corrMode', 'pos');
% STOCSY2D(H, X, ppmH, ppmX, pcaH.loadings(1,:), pcaX.loadings(1,:), 'numPointsH', 2048, 'numPointsX', 2048);
% *********************
%
% Written by Chaevien Clendinen & Gregory Stupp
% Edison Lab, University of Florida
% Department of Biochemistry and Molecular Biology
%
%% Parse Arguments
p = inputParser;
addRequired(p,'H',@(x)validateattributes(x,{'numeric'}, {'2d'}));
addRequired(p,'X',@(x)validateattributes(x,{'numeric'}, {'2d'}));
addRequired(p,'ppmH',@(x)validateattributes(x,{'numeric'}, {'row'}));
addRequired(p,'ppmX',@(x)validateattributes(x,{'numeric'}, {'row'}));

addOptional(p, 'loadingsH', [], @(x)validateattributes(x,{'numeric'}, {'row'}));
addOptional(p, 'loadingsX', [], @(x)validateattributes(x,{'numeric'}, {'row'}));
if exist('loadingsH','var')~=1
    loadingsH=0;
end
if exist('loadingsX','var')~=1
    loadingsX=0;
end

addParamValue(p, 'overlayHSQC', '');
addParamValue(p, 'threshH', 20, @(x)validateattributes(x,{'numeric'}, {'scalar'}));
addParamValue(p, 'threshX', 8e5, @(x)validateattributes(x,{'numeric'}, {'scalar'}));
addParamValue(p, 'numPointsH', 8192, @(x)validateattributes(x,{'numeric'}, {'scalar'}));
addParamValue(p, 'numPointsX', 16384, @(x)validateattributes(x,{'numeric'}, {'scalar'}));
addParamValue(p, 'threshCorr', 0, @(x)validateattributes(x,{'numeric'}, {'scalar'}));
addParamValue(p, 'corrMode', 'both', @(x) any(validatestring(x,{'pos','neg','both'})));
addParamValue(p, 'peakPickX', 1, @(x) (x==0 | x==1));
addParamValue(p, 'peakPickH', 1, @(x) (x==0 | x==1));
addParamValue(p, 'showplot', 0, @(x) (x==0 | x==1));

parse(p,H,X,ppmH,ppmX,loadingsH,loadingsX, varargin{:})

assert(size(H,2)==length(ppmH),'H must have rows as samples, columns as features (ppm)\nH and ppmH must be same length',1)
assert(size(X,2)==length(ppmX),'X must have rows as samples, columns as features (ppm)\nX and ppmX must be same length',1)

cellfun(@(f) evalin('caller',[f ' = p.Results.' f ';']), fieldnames(p.Results))
corrMode = validatestring(corrMode,{'pos','neg','both'});


%% Save orig hr data for plotting on side
ppmH_fr=ppmH;
H_fr=H;
ppmX_fr=ppmX;
X_fr=X;

%% Peak pick data if requested or lower resolution
if peakPickH
    [H,~,I_b,~]=opt_bucket(ppmH_fr,H_fr,.01,.2);
    ppmH = fliplr(mean(I_b,2)');
    H_rr=zeros(size(H,1),numPointsH);
    if size(H_fr,2)>numPointsH
        [ppmH_rr, ~] = downSample(ppmH_fr,H_fr,numPointsH);
    end
    HtoWZidx=[];
    for i=1:length(ppmH)
        [~,HtoWZidx(i)]=min(abs(ppmH(i)-ppmH_rr));
        H_rr(:,HtoWZidx(i))=H(:,i);
    end
else %or Lower resolution to numPointsH
    if size(H,2)>numPointsH
        [ppmH, H] = downSample(ppmH,H,numPointsH);
    end
    % Threshold data
    % remove columns (set them all == 0) where everything in the column is below thresh
    H(:,sum(H<threshH)==size(H,1))=0;
    
    %only calculate corr and covar on nonzeros columns
    H_rr=H; %H with zeros
    ppmH_rr = ppmH;
    HtoWZidx = logical(sum(logical(H)));
    H = H(:,HtoWZidx);
    ppmH = ppmH(HtoWZidx);
    HtoWZidx=find(HtoWZidx);
end
if peakPickX
    X=peakPickCarbon(ppmX_fr, X_fr, 'std');
    idx=~sum(X);
    X(:,idx)=[];
    ppmX=ppmX_fr;
    ppmX(idx)=[];
    X_rr=zeros(size(X,1),numPointsX);
    if size(X_fr,2)>numPointsX
        [ppmX_rr, ~] = downSample(ppmX_fr,X_fr,numPointsX);
    end
    XtoWZidx=[];
    for i=1:length(ppmX)
        [~,XtoWZidx(i)]=min(abs(ppmX(i)-ppmX_rr));
        X_rr(:,XtoWZidx(i))=X(:,i);
    end
else % Lower resolution to numPointsX
    if size(X,2)>numPointsX
        [ppmX, X] = downSample(ppmX,X,numPointsX);
    end
    X(:,sum(X<threshX)==size(X,1))=0;
    X_rr=X;
    ppmX_rr = ppmX;
    XtoWZidx = logical(sum(logical(X)));
    X = X(:,XtoWZidx);
    ppmX = ppmX(XtoWZidx);
    XtoWZidx=find(XtoWZidx);
end

%% Build Correlation and Covariance Matrix
tic
corrXH=(zscore(H)'*zscore(X))./(size(X,1)-1);
corrHX=(zscore(X)'*zscore(H))./(size(H,1)-1);
corrXH_fr=(zscore(H_fr)'*zscore(X))./(size(X,1)-1);
 corrHX_fr=(zscore(X_fr)'*zscore(H))./(size(H,1)-1);
% corrH=(zscore(H)'*zscore(H))./(size(H,1)-1);
% corrH_fr=(zscore(H_fr)'*zscore(H))./(size(H,1)-1);
corrX=(zscore(X)'*zscore(X))./(size(X,1)-1);
corrX_fr=(zscore(X_fr)'*zscore(X))./(size(X,1)-1);
toc

covarXH=bsxfun(@minus,H,mean(H))'*bsxfun(@minus,X,mean(X))./(size(zscore(X),1)-1);
covarHX=bsxfun(@minus,X,mean(X))'*bsxfun(@minus,H,mean(H))./(size(zscore(H),1)-1);
covarXH_fr=bsxfun(@minus,H_fr,mean(H_fr))'*bsxfun(@minus,X,mean(X))./(size(zscore(X),1)-1);
 covarHX_fr=bsxfun(@minus,X_fr,mean(X_fr))'*bsxfun(@minus,H,mean(H))./(size(zscore(H),1)-1);
% covarH=bsxfun(@minus,H,mean(H))'*bsxfun(@minus,H,mean(H))./(size(zscore(H),1)-1);
% covarH_fr=bsxfun(@minus,H_fr,mean(H_fr))'*bsxfun(@minus,H,mean(H))./(size(zscore(H),1)-1);
covarX=bsxfun(@minus,X,mean(X))'*bsxfun(@minus,X,mean(X))./(size(zscore(X),1)-1);
covarX_fr=bsxfun(@minus,X_fr,mean(X_fr))'*bsxfun(@minus,X,mean(X))./(size(zscore(X),1)-1);
toc

% %% Replace the values into the reduced res, save into output structure
% % if peakPickX
    S.corrXpp=corrX;
    S.Xpp=X;
    S.covarXpp=covarX;
    S.ppmXpp=ppmX;
% else
%     S.Xrr=X_rr;
%     S.corrXrr=corrX;
%     S.covarXrr=covarX;
%     S.ppmXrr=ppmX_rr;
% end
% 
% if peakPickH
%     S.corrHpp=corrH;
    S.ppmHpp=ppmH;
    S.Hpp=H;
%     S.covarHpp=covarH;
% else
%     S.covarHrr=covarH;
%     S.corrHrr=corrH;
%     S.Hrr=H_rr;
%     S.ppmHrr=ppmH_rr;
% end

% if peakPickH || peakpickX
    S.corrHXpp=corrHX;
    S.covarHXpp=covarHX;
    S.corrXHpp=corrXH;
    S.covarXHpp=covarXH;
% else
%     S.corrHXrr=corrHX;
%     S.covarHXrr=covarHX;
% end

% Full Resolution Data
S.Xfr=X_fr;
S.Hfr=H_fr;
S.ppmfr=ppmX_fr;
S.ppmHfr=ppmH_fr;
S.corrHXfr=corrHX_fr;
 S.corrXHfr=corrXH_fr;
% S.corrHfr=corrH_fr;
S.corrXfr=corrX_fr;
S.covarHXfr=covarHX_fr;
 S.covarXHfr=covarXH_fr;
% S.covarHfr=covarH_fr;
S.covarXfr=covarX_fr;

%for ease of plotting
% P.A=zeros(numPointsH,numPointsX);
% P.B=zeros(numPointsX);
% P.C=zeros(numPointsH);
% P.a=zeros(numPointsH,numPointsX);
% P.b=zeros(numPointsX);
% P.c=zeros(numPointsH);
% P.A(HtoWZidx,XtoWZidx)=corrHX;
% P.B(XtoWZidx,XtoWZidx)=corrX;
% P.C(HtoWZidx,HtoWZidx)=corrH;
% P.a(HtoWZidx,XtoWZidx)=covarHX;
% P.b(XtoWZidx,XtoWZidx)=covarX;
% P.c(HtoWZidx,HtoWZidx)=covarH;
toc

if showplot==1
    
    % Plot the loadings
    if length(loadingsH)>1 && length(loadingsX)>1
        
        %Downsample the loadings
        if length(loadingsH)>numPointsH
            [~, loadingsH] = downSample(ppmH_fr,loadingsH,numPointsH);
        end
        if length(loadingsX)>numPointsX
            [~, loadingsX] = downSample(ppmX_fr,loadingsX,numPointsX);
        end
        
        cmap=jet(100);
        
        %X
        corrX=loadingsX;
        covarX=loadingsX.*std(X);
        rangeX=[-1*max(abs(corrX)),max(abs(corrX))];
        linesX=NaN(size(cmap,1),size(corrX,2));
        ind=1;
        for k=rangeX(1):(rangeX(2)-rangeX(1))/size(cmap,1):rangeX(2)
            linesX(ind,(corrX>k))=covarX((corrX>k));
            ind=ind+1;
        end
        
        %H
        corrH=loadingsH;
        covarH=loadingsH.*std(H);
        rangeH=[-1*max(abs(corrH)),max(abs(corrH))];
        linesH=NaN(size(cmap,1),size(corrH,2));
        ind=1;
        for k=rangeH(1):(rangeH(2)-rangeH(1))/size(cmap,1):rangeH(2)
            linesH(ind,(corrH>k))=covarH((corrH>k));
            ind=ind+1;
        end
        
        %% Plot 
       
            %X
            [plotppmX,~]=makePrettyLinePlot(ppmX,linesX(1,:));
            for i=1:size(linesX,1)
                [~,linesX2(i,:)]=makePrettyLinePlot(ppmX,linesX(i,:));
            end
            linesX=linesX2;
            %H
            [plotppmH,~]=makePrettyLinePlot(ppmH,linesH(1,:));
            for i=1:size(linesH,1)
                [~,linesH2(i,:)]=makePrettyLinePlot(ppmH,linesH(i,:));
            end
            linesH=linesH2;
       
        
        % Plotting X on Y axis
        figure;
        h.ax(1)=subplot(5,5,[1 6 11 16]);
        plot(linesX(1,:),plotppmX,'Color',cmap(1,:))
        hold on
        for k=2:size(cmap,1)
            plot(linesX(k,:),plotppmX,'Color',cmap(k,:));
        end
        set(gca,'YDir','rev');
        xlim([min(min(linesX)), max(max(linesX))])
        ylim([0, max(plotppmX)])
        ylabel('Chemical Shift')
        xlabel('Correlations/Loadings')
        caxis([rangeX(1) rangeX(2)])
        
        % Plotting H on X axis
        h.ax(2)=subplot(5,5,[22 23 24 25]);
        plot(plotppmH,linesH(1,:),'Color',cmap(1,:))
        hold on
        for k=2:size(cmap,1)
            plot(plotppmH,linesH(k,:),'Color',cmap(k,:));
        end
        set(gca,'XDir','rev');
        ylim([min(min(linesH)), max(max((linesH)))])
        xlim([min(plotppmH) max(plotppmH)])
        xlabel('Chemical Shift')
        ylabel('Correlations/Loadings')
        caxis([rangeH(1) rangeH(2)])
        
    else %don't plot loadings. Plot spectra
        %X
        figure;
        h.ax(1)=subplot(5,5,[1 6 11 16]);
        plot(X_rr,ppmX_rr)
        set(gca,'YDir','rev');
        xlim([min(min(X_rr)), max(max(X_rr))])
        ylim([0, max(ppmX_rr)])
        ylabel('Chemical Shift')
        xlabel('Correlations/Loadings')
        %H
        h.ax(2)=subplot(5,5,[22 23 24 25]);
        plot(ppmH_fr,H_fr)
        set(gca,'XDir','rev');
        ylim([min(min(H_fr)), max(max(H_fr))])
        xlim([0, max(ppmH_fr)])
        xlabel('Chemical Shift')
        ylabel('Correlations/Loadings')
    end
    
    %%
    corrHXplot=P.A; %S.corrHXrr;
    corrHXplot(abs(corrHXplot)<threshCorr)=0;
    if strcmp(corrMode,'pos')
        corrHXplot(corrHXplot<0)=0;
    elseif strcmp(corrMode,'neg')
        corrHXplot(corrHXplot>0)=0;
    end
    
    
    %%
    if ~isempty(overlayHSQC)
        p=perimHSQC(overlayHSQC);
        hsqc = overlayImagesc(p.ppmH,ppmH,p.ppmX,ppmX,p.perim,corrHXplot);
        hsqc.z2=hsqc.z2>.1;
        
        corrRGB=ind2rgb(gray2ind(mat2gray(hsqc.z1),255),redblue(255));
        overlay=imoverlay(corrRGB,hsqc.z2,[0 0 0]);
        h.ax(3)=subplot(5,5,[2:5 7:10 12:15 17:20]);
        imagesc(hsqc.ppmH,hsqc.ppmX,overlay,'Parent',h.ax(3))
        set(gca,'XDir','rev');
        ylim([0 max(hsqc.ppmX(1,:))])
        xlim([min(hsqc.ppmH) max(hsqc.ppmH)])
        
        caxis([min(min(corrHXplot)) max(max(corrHXplot))])
    else
        h.ax(3)=subplot(5,5,[2:5 7:10 12:15 17:20]);
        imagesc(ppmH_rr,ppmX_rr,corrHXplot','Parent',h.ax(3))
        
        set(gca,'XDir','rev');
        ylim([0 max(ppmX_rr(1,:))])
        xlim([min(ppmH_rr) max(ppmH_rr)])
        
        caxis([min(min(corrHXplot)) max(max(corrHXplot))])
    end
    
    if strcmp(corrMode,'pos')
        colormap(blue)
    elseif strcmp(corrMode,'neg')
        colormap(red)
    else
        colormap(redblue)
    end
    
    h.link(1) = linkprop([h.ax(1),h.ax(3)],'YLim');
    h.link(2) = linkprop([h.ax(2) h.ax(3)],'XLim');
    assignin('base','h',h)
elseif showplot==0
end
end
 
function [ppm2,int2]=makePrettyLinePlot(ppm,int)
ppm2=[];
int2=[];
for i=1:length(ppm)
    ppm2=[ppm2,ppm(i)-0.01,ppm(i),ppm(i)+0.01];
    int2=[int2,0,int(i),0];
end
end
 
function [ppmL, XL] = downSample(ppm,X,numPoints)
 
ppmL = linspace(min(ppm),max(ppm),numPoints);
for i=1:size(X,1)
    XL(i,:) = interp1(ppm, X(i,:), ppmL);
end
end
