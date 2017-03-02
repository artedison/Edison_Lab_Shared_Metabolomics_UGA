% QA/QC script

%% get data ready!!!!
%% load in spectra

loadallft
%% remove the redos unless using the fixed list

% spectra(11) = [];
% spectra(20) = [];
% spectra(21) = [];
% spectra(23) = [];
% spectra(47) = [];
%  spectra(57) = [];
%  spectra(63) = [];
% spectra(71) = [];

%% reference the spectra

[spectra] = ref_spectra(spectra,0.01)


%% create the sample matrix

[X,ppm,XTitles]=Setup1D(spectra);

%%

Xtitles=XTitles

%% removes anything but numbers and changes string to number

for i = 1:length(XTitles)
    s = XTitles{i};
    XTitles(i) = {s(1:end-3)};
end

XTitles=cellfun(@str2num, XTitles)

%% visulize the spectra

plotr(ppm,X)

%%
whichLine()

%%
highlightLine(XTitles)

%% flag just look wierd as a spectroscopist

questionable=[64]

flag_different_by_eye=ismember(XTitles,questionable)

%% Yvec=  making the y (identity vector) 
% 0 = ct/5 2=control 13=natl history 14= pooled 18=steroid control
% 19=steroid treatment 15=blanks

% need to fix or make sure that the groups do not start with same letter
% otherwise is made into the same group

[Yvec]= group2('BR1_BLM_2_Core2_sheet_2.xlsx')

%% remove water and ends

XR=remove_region(X,ppm,4.8,5.15);
[XR1,ppm]=remove_ends(XR,ppm,-0.2,10.0);


%% checked the difference in the temperatures
figure
hold on
plotr(ppm,XR1)
hold off

%% align spectra


XAL=guide_align1D(XR1,ppm,'correlation','CCOW'); % USE!! takes a while but
% still the best
% XAL2=guide_align1D(XR1,ppm,'correlation','PAFFT');
%XAL3=guide_align1D(XR1,ppm,'correlation','RAFFT');

%% and visulaize

figure 
hold on
plotr(ppm,XAL);
hold off

%%
whichLine()

%%
highlightLine(XTitles)


%%

XAL_norm=XAL;

%%

XAL_norm(70,:)=[]; %blank
% XAL_norm(63,:)=[]; %RG/dilution issue
% XAL_norm(57,:)=[]; %RG/dilution issue
% XAL_norm(46,:)=[]; %blank
XAL_norm(29,:)=[]; %blank


%%

Yvec1=Yvec

%%
Yvec(70)=[];
% Yvec(63)=[];
% Yvec(57)=[];
% Yvec(46)=[];
Yvec(29)=[];

%%

XTitles1=XTitles
%%

Xtitles1=Xtitles

%%
XTitles(70)=[];
XTitles(29)=[];
% XTitles(57)=[];
% XTitles(46)=[];
% XTitles(1)=[];

%%

Xtitles(70)=[];
Xtitles(29)=[];

%% Normalize??... this is controversial.. in theory?


%%

XALN=normalize(XAL_norm,ppm,'PQN');

normcheck(XALN)

%%

figure 
hold on 
plotr(ppm,XALN)
hold off

%%

whichLine()

%%

highlightLine(Xtitles)

%% create a flag vector of those that are questionable data
% where 1=wonky and 0=normal

questionable_norm=[6,18,22,31,33,29,43,50,54,69];

flag_norm=ismember(XTitles,questionable_DSS_int);


%% 0 = ct/5 2=control 13=natl history 14= pooled 18=steroid control
% 19=steroid treatment 15= blanks

%X_blanks=XALN(Yvec==15,:); 
X_std_treat=XALN(Yvec==19,:);
X_std_con=XALN(Yvec==18,:);
X_pooled=XALN(Yvec==14,:);
X_natl=XALN(Yvec==13,:);
X_con=XALN(Yvec==2,:);
X_ct=XALN(Yvec==0,:);
X_dmd=vertcat(X_std_con,X_std_treat,X_natl,X_ct);

%% 0 = ct/5 2=control 13=natl history 14= pooled 18=steroid control
% 19=steroid treatment 15= blanks

%X_blanks_Titles=XTitles(Yvec==15,:); 
X_std_treat_Titles=XTitles(Yvec==19,:);
X_std_con_Titles=XTitles(Yvec==18,:);
X_pooled_Titles=XTitles(Yvec==14,:);
X_natl_Titles=XTitles(Yvec==13,:);
X_con_Titles=XTitles(Yvec==2,:);
X_ct_Titles=XTitles(Yvec==0,:);
X_dmd_Titles=vertcat(X_std_con_Titles,X_std_treat_Titles,X_natl_Titles,X_ct_Titles);

%% 0 = ct/5 2=control 13=natl history 14= pooled 18=steroid control
% 19=steroid treatment 15= blanks

%X_blanks_titles=Xtitles(Yvec==15,:); 
X_std_treat_titles=Xtitles(Yvec==19,:);
X_std_con_titles=Xtitles(Yvec==18,:);
X_pooled_titles=Xtitles(Yvec==14,:);
X_natl_titles=Xtitles(Yvec==13,:);
X_con_titles=Xtitles(Yvec==2,:);
X_ct_titles=Xtitles(Yvec==0,:);
X_dmd_titles=vertcat(X_std_con_titles,X_std_treat_titles,X_natl_titles,X_ct_titles);


%% remove the wonky

% XALN(90,:)=[];
% XALN(63,:)=[];
% XALN(57,:)=[];
% XALN(46,:)=[];
% XALN(1,:)=[];



%% Start QC 

%% A nice overview of your data. Displays histogram and box plots of log-fold change vs. median for all
% features in each spectrum.  Dilution / normalization effects are often
% visible as distributions not centered at 0.

normcheck(XAL)

%% create a flag vector of those that are questionable data
% where 1=wonky and 0=normal

questionable_variance=[29,70]
flag_variance=ismember(XTitles,questionable_variance)
 
%% checked that temperature hasn't changed... load in xls from greg on terminal

figure
hold on
scatter(temperature(:,1),temperature(:,2))
text(temperature(:,1),temperature(:,2),Xtitles)
hold off


%% create a flag vector of those that are questionable data
% where 1=wonky and 0=normal
% wrong vector!!! REDO these are indeces
questionable_temperature=[56,58,59,60,61,62,64];

flag_temperature=ismember(XTitles,questionable_temperature);

%% look at the variance in temperature effects spectra before alignment
figure
hold on
plotr(ppm,XR1(1:20,:)+10000000,'b')
plotr(ppm,XR1(20:end,:),'r')
hold off

%% look at the variance in temperature effects spectra after alignment
figure
hold on
plotr(ppm,XAL(1:20,:)+10000000,'b')
plotr(ppm,XAL(20:end,:),'r')
hold off

%% integrate DSS before normalization

DSS_integration=IntegralPeak_roi(XAL,ppm,-0.04,0.04);

%%
figure
hold on
scatter(Yvec,DSS_integration)
text(Yvec+1,DSS_integration,Xtitles,'FontSize', 15)
title('DSS integration distribution per group')
    xlabel('Group')
    ylabel('Integration')
hold off

%% create a flag vector of those that are questionable data
% where 1=wonky and 0=normal

questionable_DSS_int=[89,35,39,75,71,52,50,23,64];

flag_DSS_int=ismember(XTitles,questionable_DSS_int)


%% normalize?? difference? how to choose which one.. decide if what you are seeing is because of
% instrument variation or not.. technical ok.. others are more concerning

DSS_integration_postnorm=IntegralPeak_roi(XALN,ppm,-0.04,0.04);


%% separate into groups 0 = ct/5 2=control 13=natl history 15= pooled 18=steroid control
% 19=steroid treatment 1= blanks

figure
hold on
scatter(Yvec,DSS_integration_postnorm)
text(Yvec,DSS_integration_postnorm,Xtitles,'FontSize', 15)
title('Normalized DSS integration distribution per group')
    xlabel('Group')
    ylabel('Integration')
hold off

%% create a flag vector of those that are questionable data
% where 1=wonky and 0=normal

questionable_nDSS_int=[89,71,35,50,67,39,75,52,50,23,73];

flag_nDSS_int=ismember(XTitles,questionable_nDSS_int);


%% FWHM DSS before norm

[Fullw_DSS] = fwhm(XAL,ppm,-0.04,0.04)

%%

figure
hold on
scatter(Yvec1,Fullw_DSS)
text(Yvec1+1,Fullw_DSS,Xtitles1,'FontSize',15)
title('Full Width of the Half Max distribution on DSS per group')
    xlabel('Group')
    ylabel('FWHM')
hold off

%% flag or decide first?? if biological or instrumentation?

%% FWHM DSS before norm

questionable_DSS_fwhm=[29,66,86,45];

flag_DSS_fwhm=ismember(XTitles,questionable_DSS_fwhm);

%% FWHM DSS before norm

[Fullw_DSS_norm] = fwhm(XALN,ppm,-0.04,0.04)

%%

figure
hold on
scatter(Yvec,Fullw_DSS_norm)
text(Yvec,Fullw_DSS_norm,Xtitles,'FontSize',15)
title('Full Width of the Half Max distribution on normalized DSS per group')
    xlabel('Group')
    ylabel('FWHM')
hold off


%% create a flag vector of those that are questionable data
% where 1=wonky and 0=normal

questionable_nDSS_fwhm=[29,66,86,45];

flag_nDSS_fwhm=ismember(XTitles,questionable_nDSS_fwhm);

%% peakpicking for BA plots?

[pooled_peaks,pooled_shifts]=Peakpick1D(X_pooled,ppm,'mean',0.25);
[con_peaks,con_shifts]=Peakpick1D(X_con,ppm,'mean',0.25);
[dmd_peaks,dmd_shifts]=Peakpick1D(X_dmd,ppm,'mean',0.25);

%%

idx=[1:length(A)]
ind=num2cell(idx)
idx=cellfun(@num2str, ind,'UniformOutput' , false)

%%

[peaks,shifts]=Peakpick1D(XALN,ppm,'mean',0.2);


%%

[cr, fig, statsStruct] = BlandAltman(pooled_peaks(1,:)',pooled_peaks(4,:)');
%BlandAltman(data1, data2,label,tit,gnames,corrinfo,BAinfo,limits,colors,symbols)
%text(pooled_shifts,pooled_peaks)

%%

whichLine()

%%

highlightLine(idx)
%%

[cr, fig, statsStruct] = BlandAltman(X_pooled(1,:)',X_pooled(5,:)');


%% blanks

autoBA(X_blanks,X_blanks_titles)

%% create a flag vector of those that are questionable data
% where 1=wonky and 0=normal

questionable_BA_blanks=[46,90];

flag_BA_blanks=ismember(XTitles,questionable_BA_blanks);

%% pooled

autoBA(pooled_peaks,X_pooled_titles)

%% create a flag vector of those that are questionable data
% where 1=wonky and 0=normal

questionable_DSS_fwhm=[fine];

flag_DSS_fwhm=ismember(XTitles,questionable_DSS_fwhm);

%% controls

autoBA(con_peaks,X_con_titles)

%% create a flag vector of those that are questionable data
% where 1=wonky and 0=normal

questionable_DSS_fwhm=[?];

flag_DSS_fwhm=ismember(XTitles,questionable_DSS_fwhm);

%% DMD

autoBA(X_dmd)

%% create a flag vector of those that are questionable data
% where 1=wonky and 0=normal

questionable_DSS_fwhm=[?];

flag_DSS_fwhm=ismember(XTitles,questionable_DSS_fwhm);

%% and visulaize
figure 
hold on
plotr(ppm,XALN);
hold off

%%
whichLine()

%%
highlightLine(XTitles)

%% integrate creatine

creatine=IntegralPeak_roi(XALN,ppm,3.92,3.94);

%% integrate creatinine

creatinine=IntegralPeak_roi(XALN,ppm,4.04,4.065);

%% create ratio

cc_ratio=creatine./creatinine

%% plot distribution  0 = ct/5 2=control 13=natl history 15= pooled 18=steroid control
% 19=steroid treatment 1= blanks

figure
hold on
scatter(Yvec,cc_ratio)
text(Yvec+0.5,cc_ratio,Xtitles,'FontSize',15)
title('Distibution of creatine/creatinine ratio')
    xlabel('Group')
    ylabel('creatine/creatinine ratio')
hold off

%%
F=cc_ratio./repmat(median(cc_ratio),[size(cc_ratio,1),1]);

%%
figure, boxplot(F,'plotstyle','compact','symbol',' ')
    %ylim([-4 4])
    title('Boxplots of log-fold-change vs. median values')
    xlabel('Observation')
    ylabel('log-fold-change vs. median  feature value')

%% flag problems

questionable_cc_ratio=[1,46,90];

flag_cc_ratio=ismember(XTitles,questionable_cc_ratio);


