

%% Calculate p value of correlation
% Two tailed significance testing
rho=corr;
rhocovar=covar;
n=size(H,1);

%PVAL of PEARSON Tail probability for Pearson's linear correlation.
t = rho.*sqrt((n-2)./(1-rho.^2)); % +/- Inf where rho == 1
pval = 2*tcdf(-abs(t),n-2);

%% Get the indeces for all the insignificant correlations and remove them from the correlation matrix

idx=find(pval>0.05);

rho(idx)=0;
rhocovar(idx)=0;
%% Remove negative correlations
idx=find(rho<0);

rho(idx)=0;
rhocovar(idx)=0;

%idx=find(rho~=0);
%% find the peaks that are not zero (positively significant)
[r,c]=find(rho>0);

%% Remove repeating chem shifts

c=unique(c);
r=unique(r);


%% Just for visualization, serves no particular function in the rest of the script

ppmsigH=ppmHs(r);

Hsig=Hs(:,r);

ppmsigX=ppmXs(c);

Xsig=Xs(:,c);

%% Put covariance on a 0 to 1 scale by dividing by the highest positive covariace overall (optional)
%rhocovar=rhocovar./max(max(rhocovar));

%% Create a PairsX2H where each structure is a carbon chem shift annd contains the chem shifts of protons highly correlated
for i=1:length(c)
PairsX2H(i).ppm_of_carbon=ppmXs(c(i)); %the carbon chem shift

PairsX2H(i).idxH=find(rho(:,c(i))~=0); %indeces of the proton chem shifts

PairsX2H(i).covarH=rhocovar(PairsX2H(i).idxH,c(i))./max(max(rhocovar(PairsX2H(i).idxH,c(i)))); %covariance of the carbon to each H shift

PairsX2H(i).ppmH=ppmHs(PairsX2H(i).idxH); % the proton chem shifts
end

%% Pick out peaks that are greater than half of the maximum covariance and put them into a matrix called ppmcovar
for i=1:length(PairsX2H)
PairsX2H(i).topcovar=find((ceil(PairsX2H(i).covarH * 100)./100)>0.5); %the indeces of the top peaks
PairsX2H(i).ppmcovar=PairsX2H(i).ppmH(PairsX2H(i).topcovar);
end

%% Combining ppm shifts within 0.3 range

for i=1:length(PairsX2H)
PairsX2H(i).ppmcovar=comchem(PairsX2H(i),'ppmcovar',0.03);
end

%% Find C shifts that share one or more ppmcovar H shifts
for i=1:length(PairsX2H)
    for j=1:length(PairsX2H(i).ppmcovar)
        for p=1:length(PairsX2H)
            for d=1:length(PairsX2H(p).ppmcovar)
                if ~isempty(find(PairsX2H(i).ppmcovar(j)>=PairsX2H(p).ppmcovar(d)-0.03)) && ~isempty(find(PairsX2H(i).ppmcovar(j)<=PairsX2H(p).ppmcovar(d)+.03))
                    PairsX2H(i).cshifts(p)=PairsX2H(p).ppm_of_carbon;
                else
                    PairsX2H(i).cshifts(p)=0;
                end
            end
        end
    end
end

for i=1:length(PairsX2H);
    PairsX2H(i).idC=find(PairsX2H(i).cshifts~=0);
    PairsX2H(i).cshifts(PairsX2H(i).cshifts==0)=[];
end
    
%% Combine all common C and H shifts into a structure

for i=1:length(PairsX2H)
    for j=1:length(PairsX2H(i).idC)
        Plist(i).C(j)=PairsX2H(PairsX2H(i).idC(j)).ppm_of_carbon;
    end
end

p=1;
for i=1:length(PairsX2H)
    for j=1:length(PairsX2H(i).idC)
        Plist(i).H(p:p+length((PairsX2H(PairsX2H(i).idC(j)).ppmcovar))-1)=PairsX2H(PairsX2H(i).idC(j)).ppmcovar;

        p=length(Plist(i).H)+1;
    end
    Plist(i).H=unique(Plist(i).H);
    p=1;
end

        
%% combine similar chemical shifts for the H
for i=1:length(Plist)
Plist(i).H=comchem(Plist(i),'H',0.03);
end

%% Remove duplicate peaklist

for i=1:length(PairsX2H)
    for j=2:length(PairsX2H(i).idC)
        if isequal(Plist(i),Plist(PairsX2H(i).idC(j)))
            Plist(PairsX2H(i).idC(j))=setfield(Plist(PairsX2H(i).idC(j)),'C',NaN);
        end
    end
end

i=1;
while i<length(Plist)
    if isnan(Plist(i).C)
        Plist(i)=[];
    else
        i=i+1;
    end
end

%%
