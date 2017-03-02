function X = peakPickCarbon(ppmX, X, thresh, showplot)

% Thresh is either a number (absolute value), 'std', or 'iqr'

%%
if exist('showplot','var')~=1
    showplot=0;
end

%%
ppmError=0.3;
for spec=1:size(X,1)
    y=X(spec,:);
    
    %determine thresh
    if ~isnumeric(thresh)
        if strcmpi(thresh,'std')
            t = std(y)*2;
        elseif strcmpi(thresh, 'iqr')
            t = abs(quantile(y,.75)-3.0*(quantile(y,.75)-quantile(y,.25)));
        else
            error('Thresh is either a number (absolute value), ''std'', or ''iqr''')
        end
    elseif length(thresh)>1
        t = thresh(spec);
    else
        t = thresh;
    end
        
    y(y<t)=0;
    if sum(y)==0
        error(['Your thresh is too high. Sample: ' num2str(spec)])
    end
    [maxtab, ~] = peakdet(y,1);
    % merge together peaks within 0.2 ppm
    midx=find(diff(ppmX(maxtab(:,1)))<=0.2);
    for i=1:length(midx)
        maxtab(midx(i),2)=sum(maxtab(midx(i):midx(i)+1,2));
        maxtab(midx(i),1)=round(mean((maxtab(midx(i):midx(i)+1,1))));
        maxtab(midx(i)+1,:)=[];
        midx=midx-1;
    end
    peaks{spec}=maxtab(:,1);
    ppm{spec}=ppmX(maxtab(:,1));
end
d=ppmError/mean(diff(ppmX));
[peakGroup,ppmIdx,sampleIdx]=groupPeaks(peaks, d);

X_o=X;
X=zeros(size(X));
%%
for i=1:max(peakGroup)
    p=ppmIdx(find(peakGroup==i));
    s=sampleIdx(find(peakGroup==i));
    [~,idx]=min(abs(ppmX-ppmX(round(mean(p)))));
    p2=[];
    p2(s)=p;
    p=p2;
    for sample=s
        X(sample,idx)=X_o(sample,p(sample));
    end
end
%%
if showplot
figure,plot(ppmX,X), h1=gca;
figure,plot(ppmX,X_o), h2=gca;
linkaxes([h1, h2])
end