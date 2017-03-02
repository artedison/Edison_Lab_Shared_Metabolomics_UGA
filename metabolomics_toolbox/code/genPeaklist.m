function [Plist,PairsX2H]=genPeaklist(corr,covar,replicates,ppmHs,ppmXs,varargin)
%%*****************

%% Parse Arguments
m = inputParser;
addRequired(m,'corr',@(x)validateattributes(x,{'numeric'}, {'2d'}));
addRequired(m,'covar',@(x)validateattributes(x,{'numeric'}, {'2d'}));
addRequired(m,'replicates',@(x)validateattributes(x,{'numeric'}, {'2d'}));
addRequired(m,'ppmHs',@(x)validateattributes(x,{'numeric'}, {'row'}));
addRequired(m,'ppmXs',@(x)validateattributes(x,{'numeric'}, {'row'}));

addParamValue(m, 'PvalThreshold',0.05, @(x)validateattributes(x,{'numeric'}, {'scalar'}));
addParamValue(m, 'FDR',1, @(x) (x==0 | x==1 | x==2));
addParamValue(m, 'CovarThreshold', 0.7, @(x)validateattributes(x,{'numeric'}, {'scalar'}));
addParamValue(m, 'ShiftRange', [0.05 0.05], @(x)validateattributes(x,{'numeric'}, {'size',[1,2]}));
addParamValue(m, 'options', 2, @(x) (x==1 | x==2));
addParamValue(m, 'common', 5, @(x)validateattributes(x,{'numeric'}, {'scalar'}));
addParamValue(m, 'homo', 1, @(x) (x==1 | x==0));


%option 1 is for corr first then covar 2

parse(m,corr,covar,replicates,ppmHs,ppmXs,varargin{:})
cellfun(@(f) evalin('caller',[f ' = m.Results.' f ';']), fieldnames(m.Results))


%% Remove negative correlations

idx=find(corr<0);

corr(idx)=0;
covar(idx)=0;

%idx=find(rho~=0);
%% find the peaks that are not zero (positively significant)
[~,c]=find(corr>0);

%% Remove repeating chem shifts

c=unique(c);

%% Put covariance on a 0 to 1 scale by dividing by the highest positive covariace overall (optional)
% %rhocovar=rhocovar./max(max(rhocovar));

%% Create a PairsX2H where each structure is a carbon chem shift annd contains the chem shifts of protons highly correlated

for i=1:length(c)
    PairsX2H(i).ppm_of_carbon=ppmXs(c(i)); %the x chem shift
    
    PairsX2H(i).idxH=find(corr(:,c(i))>0); %indeces of the proton chem shifts
    
    PairsX2H(i).covarH=covar(PairsX2H(i).idxH,c(i))./max(max(covar(PairsX2H(i).idxH,c(i)))); %covariance of the carbon to each H shift
    %     PairsX2H(i).covarH=covar(PairsX2H(i).idxH,c(i)); %covariance of the carbon to each H shift
    
    PairsX2H(i).corrH=corr(PairsX2H(i).idxH,c(i))./max(max(corr(PairsX2H(i).idxH,c(i)))); %correlation of the carbon to each H shift
    
    
end

%% Combine like Peaks across the structure
PairsX2H=com(PairsX2H,ShiftRange(2));

%% 1=corr first then covariance 2=covar first then corr
for i=1:length(PairsX2H)
    PairsX2H(i).ppmH=ppmHs(PairsX2H(i).idxH);
end

if options==1
    
    % Filter by correlation significance after fdr correction
    n=replicates;
    for i=1:length(PairsX2H)
        PairsX2H(i).t = PairsX2H(i).corrH.*sqrt((n-2)./(1-PairsX2H(i).corrH.^2));
        PairsX2H(i).pval = 2*tcdf(-abs(PairsX2H(i).t),n-2);
    end
    if FDR==1
        for i=1:length(PairsX2H)
            [~, PairsX2H(i).crit_p, PairsX2H(i).adj_p]=fdr_bh(PairsX2H(i).pval);
        end
    elseif FDR==0
        for i=1:length(PairsX2H)
            PairsX2H(i).crit_p=PvalThreshold;
        end
    elseif FDR==2
        for i=1:length(PairsX2H)
            if (size(PairsX2H(i).pval(PairsX2H(i).pval<PvalThreshold))/size(PairsX2H(i).pval))<=0.26
                PairsX2H(i).crit_p=PvalThreshold;
            else
                [~, PairsX2H(i).crit_p, PairsX2H(i).adj_p]=fdr_bh(PairsX2H(i).pval);
            end
        end
    end
    
    
    for i=1:length(PairsX2H)
        PairsX2H(i).idxpval=find(PairsX2H(i).pval<=PairsX2H(i).crit_p);
        PairsX2H(i).ppmcorr=PairsX2H(i).ppmH(PairsX2H(i).idxpval);
        PairsX2H(i).corrsH=PairsX2H(i).corrH(PairsX2H(i).idxpval);
        PairsX2H(i).covarsH=PairsX2H(i).covarH(PairsX2H(i).idxpval);
    end
    
    % divide by the highest covariance and filter based on covariances
    % greater than covarthreshold
    if homo==1 
        for i=1:length(PairsX2H)
            [~,u]=min(abs(PairsX2H(i).ppmcorr - (PairsX2H(i).ppm_of_carbon)));
            PairsX2H(i).covarthresh=[PairsX2H(i).covarsH(u)*2 PairsX2H(i).covarsH(u)/(1/CovarThreshold)];
            PairsX2H(i).topcovar=find((PairsX2H(i).covarsH)<PairsX2H(i).covarthresh(1) & (PairsX2H(i).covarsH)>PairsX2H(i).covarthresh(2));
            if isempty(PairsX2H(i).topcovar)
                continue
            end
            PairsX2H(i).toppm=PairsX2H(i).ppmcorr((PairsX2H(i).topcovar));
            PairsX2H(i).corrsH=PairsX2H(i).corrsH(PairsX2H(i).topcovar);
            PairsX2H(i).covarsH=PairsX2H(i).covarsH(PairsX2H(i).topcovar);
        end
    elseif homo==0
        for i=1:length(PairsX2H)
            PairsX2H(i).topcovar=find((ceil(PairsX2H(i).covarH * 100)./100)>CovarThreshold); %the indeces of the top peaks
            if isempty(PairsX2H(i).topcovar)
                continue
            end
            PairsX2H(i).toppm=PairsX2H(i).ppmcorr((PairsX2H(i).topcovar));
            PairsX2H(i).corrsH=PairsX2H(i).corrsH(PairsX2H(i).topcovar);
            PairsX2H(i).covarsH=PairsX2H(i).covarsH(PairsX2H(i).topcovar);
        end
    end
    
elseif options==2
    % Get covariance Threshold
    for i=1:length(PairsX2H)
        if homo==1
            [~,u]=min(abs(PairsX2H(i).ppmH - (PairsX2H(i).ppm_of_carbon)));
            PairsX2H(i).covarthresh=[PairsX2H(i).covarH(u)*2 PairsX2H(i).covarH(u)/(1/CovarThreshold)];
            PairsX2H(i).topcovar=find((PairsX2H(i).covarH)<PairsX2H(i).covarthresh(1) & (PairsX2H(i).covarH)>PairsX2H(i).covarthresh(2));
            PairsX2H(i).ppmcovar=PairsX2H(i).ppmH(PairsX2H(i).topcovar);
            PairsX2H(i).corrsH=PairsX2H(i).corrH(PairsX2H(i).topcovar);
            PairsX2H(i).covarsH=PairsX2H(i).covarH(PairsX2H(i).topcovar);
        elseif homo==0
            PairsX2H(i).topcovar=find((ceil(PairsX2H(i).covarH * 100)./100)>CovarThreshold); %the indeces of the top peaks
            PairsX2H(i).ppmcovar=PairsX2H(i).ppmH(PairsX2H(i).topcovar);
            PairsX2H(i).corrsH=PairsX2H(i).corrH(PairsX2H(i).topcovar);
            PairsX2H(i).covarsH=PairsX2H(i).covarH(PairsX2H(i).topcovar);
        end
    end
    
    % Filter by correlation significance after fdr correction
    n=replicates;
    for i=1:length(PairsX2H)
        PairsX2H(i).t = PairsX2H(i).corrsH.*sqrt((n-2)./(1-PairsX2H(i).corrsH.^2));
        PairsX2H(i).pval = 2*tcdf(-abs(PairsX2H(i).t),n-2);
    end
    if FDR==1
        for i=1:length(PairsX2H)
            [~, PairsX2H(i).crit_p, PairsX2H(i).adj_p]=fdr_bh(PairsX2H(i).pval);
        end
    elseif FDR==0
        for i=1:length(PairsX2H)
            PairsX2H(i).crit_p=PvalThreshold;
        end
    elseif FDR==2
        for i=1:length(PairsX2H)
            if (size(PairsX2H(i).pval(PairsX2H(i).pval<PvalThreshold))/size(PairsX2H(i).pval))<=0.26
                PairsX2H(i).crit_p=PvalThreshold;
            else
                [~, PairsX2H(i).crit_p, PairsX2H(i).adj_p]=fdr_bh(PairsX2H(i).pval);
            end
        end
    end
    for i=1:length(PairsX2H)
        PairsX2H(i).idxpval=find(PairsX2H(i).pval<=PairsX2H(i).crit_p);
        PairsX2H(i).idxpval1=find(PairsX2H(i).corrsH>0.5);
        si=ismember(PairsX2H(i).idxpval,PairsX2H(i).idxpval1);
        PairsX2H(i).toppm=PairsX2H(i).ppmcovar(PairsX2H(i).idxpval(si));
    end
end



%% Combining ppm shifts within specified range

for i=1:length(PairsX2H)
    PairsX2H(i).toppm=comchem(PairsX2H(i),'toppm',ShiftRange(1));
end


%% Find X shifts that share one or more ppmcovar H shifts
for i=1:length(PairsX2H)
    if isempty(PairsX2H(i).toppm)
        continue
    end
    PairsX2H(i).cshifts=[];
    PairsX2H(i).idC=[];
    
    for j=1:length(PairsX2H)
        if isempty(PairsX2H(j).toppm)
            continue
        end
        if homo==1
            for k=1:length(PairsX2H(i).toppm)
                if sum(ismember(PairsX2H(i).toppm,PairsX2H(j).toppm))>=common && any(PairsX2H(i).ppm_of_carbon-ShiftRange(1)<PairsX2H(j).toppm & PairsX2H(i).ppm_of_carbon+ShiftRange(1)>PairsX2H(j).toppm)
                    if any(PairsX2H(j).ppm_of_carbon-ShiftRange(1)<PairsX2H(i).toppm & PairsX2H(j).ppm_of_carbon+ShiftRange(1)>PairsX2H(i).toppm)
                        PairsX2H(i).cshifts=[PairsX2H(i).cshifts, PairsX2H(j).ppm_of_carbon];
                        PairsX2H(i).idC=[PairsX2H(i).idC,j];
                        break
                    end
                end
            end
        elseif homo==0
            
            if sum(ismember(PairsX2H(i).toppm,PairsX2H(j).toppm))>=common
                PairsX2H(i).cshifts=[PairsX2H(i).cshifts, PairsX2H(j).ppm_of_carbon];
                PairsX2H(i).idC=[PairsX2H(i).idC,j];
            end
        end
    end
end
for i=1:length(PairsX2H)
    PairsX2H(i).idC=unique(PairsX2H(i).idC);
end

%% Combine all common X shifts into a structure

for i=1:length(PairsX2H)
    Plist(i).X=PairsX2H(i).cshifts;
    Plist(i).idC=PairsX2H(i).idC;
end
if homo==1
    for i=1:length(Plist)
        if isempty(Plist(i).X)
            Blist(i)=1;
        end
    end
    Plist(Blist==1)=[];
    for j=1:length(Plist)
        if j>length(Plist)
            break
        end
        for q=j+1:length(Plist)
            if q>length(Plist)
                continue
            end
            if length(Plist(j).idC)==length(Plist(q).idC)
                if sum(Plist(j).idC==Plist(q).idC)==length(Plist(j).idC)
                    Plist(q)=[];
                end
            end
        end
    end
    for i=1:length(Plist)
        Plist(i).names=lookupColmar(Plist(i).X,'C');
        pause(1)
    end
    return
elseif homo==0
    
    p=1;
    for i=1:length(PairsX2H)
        if isempty(PairsX2H(i).toppm) || isempty(PairsX2H(i).cshifts)
            continue
        end
        for j=1:length(PairsX2H(i).cshifts)
            Plist(i).H(p:p+length(PairsX2H(PairsX2H(i).idC(j)).ppmcovar)-1)=PairsX2H(PairsX2H(i).idC(j)).ppmcovar;
            p=length(Plist(i).H)+1;
        end
        Plist(i).H=unique(Plist(i).H);
        p=1;
        
    end
    % Combine similar chemical shifts for the H
    for i=1:length(Plist)
        if isempty(Plist(i).H)
            continue
        end
        Plist(i).H=comchem(Plist(i),'H',ShiftRange(1));
    end
end
%
for i=1:length(Plist)
    if isempty(Plist(i).X)
        Alist(i)=1;
    end
end
Plist(Alist==1)=[];

for i=1:length(Plist)
    if i==length(Plist)
        continue
    end
    for j=i+1:length(Plist)
        if j>=length(Plist)
            break
        end
        if sum(ismember(Plist(i).H,Plist(j).H))>=common 
            Plist(i).X=[Plist(i).X,Plist(j).X];
            Plist(i).H=[Plist(i).H,Plist(j).H];
        end
        Plist(j).X=[];
    end
end

for i=1:length(Plist)
    if isempty(Plist(i).X)
        Blist(i)=1;
    end
end
if exist('Blist')
Plist(Blist==1)=[];
end
for i=1:length(Plist)
    if i>length(Plist)
        continue
    end
    Plist(i).X=comchem(Plist(i),'X',ShiftRange(1));
    Plist(i).H=comchem(Plist(i),'H',ShiftRange(1));
    
end
tic
for i=1:length(Plist)
    if length(Plist(i).X)<7
        Plist(i).namesH=lookupColmar(Plist(i).X,'H');
        pause(1)
    end
    Plist(i).namesX=lookupColmar(Plist(i).H,'C');
    pause(1)
end
toc


end


function PairsX2H=com(PairsX2H,ShiftRange)
for i=1:length(PairsX2H)
    if i==length(PairsX2H) && PairsX2H(i).ppm_of_carbon==PairsX2H(i-1).ppm_of_carbon
        continue
    elseif i==length(PairsX2H) && PairsX2H(i).ppm_of_carbon~=PairsX2H(i-1).ppm_of_carbon
        z(i)=0;
        continue
    end
    
    if PairsX2H(i).ppm_of_carbon>=PairsX2H(i+1).ppm_of_carbon-ShiftRange &&  PairsX2H(i).ppm_of_carbon<=PairsX2H(i+1).ppm_of_carbon+ShiftRange
        z(i:i+1)=[PairsX2H(i).ppm_of_carbon PairsX2H(i+1).ppm_of_carbon];
    end
end

%%
zi=1;
for j=1:length(PairsX2H)-1
    if z(j)==0
        continue
    elseif z(j)>=z(j+1)-ShiftRange && z(j)<=z(j+1)+ShiftRange
        continue
    else
        a=[z(zi:j) 0];
        id=[0 logical(a) 0];
        ii1=strfind(id,[0,1]);
        ii2=strfind(id,[1,0])-1;
        for k=1:length(ii1)
            a(ii1(k):ii2(k))=median(a(ii1(k):ii2(k)));
        end
        z(zi:j)=a(1:end-1);
        zi=j+1;
    end
end

z(zi:end)=median(z(zi:end));
%%
idx=find(z~=0);

for i=1:length(idx)
    PairsX2H(idx(i)).ppm_of_carbon=z(idx(i));
end
%%

for i=1:length(PairsX2H)
    j=1;
    if i==length(PairsX2H)
        continue
    end
    while PairsX2H(i).ppm_of_carbon==PairsX2H(i+j).ppm_of_carbon;
        
        PairsX2H(i).idxH(end+1:length(PairsX2H(i).idxH)+length(PairsX2H(i+j).idxH))=PairsX2H(i+j).idxH; %indeces of the proton chem shifts
        
        PairsX2H(i).covarH(end+1:length(PairsX2H(i).covarH)+length(PairsX2H(i+j).covarH))=PairsX2H(i+j).covarH; %covariance of the carbon to each H shift
        
        PairsX2H(i).corrH(end+1:length(PairsX2H(i).corrH)+length(PairsX2H(i+j).corrH))=PairsX2H(i+j).corrH; %correlation of the carbon to each H shift
        
        j=j+1;
    end
end

%%
for i=1:length(z)
    if z(i)==0
        z(i)=i;
    end
end

[p,iz,~]=unique(z);

iz=sort(iz);

for i=1:length(iz)
    v(i)=PairsX2H(iz(i));
end
PairsX2H=v;
end
