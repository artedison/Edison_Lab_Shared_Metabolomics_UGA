function XALg=guide_align1D(X,ppm,distance_metric,alignment_method)

    % Author: Steven L Robinette
    % Ver 0.1
    % Tested on Matlab Version R2017b
    % Date: 25FEB2019
    %
    %
    % Description:
    %   Calculates guided alignment using hierarchical clustering and an
    %   alignment method of your choice - default is PAFFT. For mention of
    %   "guide" vs. "star" alignment, see Anal Chem. 2011 Mar 1; 83(5): 1649?1657.
    %
    % Input:
    %   X: Data matrix of spectra
    %   ppm: chemical shift vector
    %   distance_metric: string, either 'correlation' or 'spearman'
    %   alignment_method: string, either 'CCOW','ICOSHIFT','RAFFT','PAFFT'
    %   NOTE: 'ICOSHIFT' does not run correctly.
    %
    % Output:
    %   XALg: Guide-aligned spectra    
    %
    % Log:
    %   Edited by : MTJ,LM,YW,SZ
    %   Date      : 25FEB2019
    %   Ver       : 0.1
    %       -renamed "clustering_method" -> "distance_metric"
    %       -ICOSHIFT method pops an error
    % Example run:
    %


XALg=X;
XALg(isnan(XALg))=0;

if exist('cluster_method')~=1
    distance_metric='correlation';
end
if exist('alignment_method')~=1
    alignment_method='PAFFT';
end

%bucket spectra in .1 ppm increments between 0-10 ppm for low-grain clustering
[h,k1]=min(abs(ppm-0));
[h,k2]=min(abs(ppm-10));
bucket_domain=round(k1:(k2-k1)/100:k2);
Xbucket=zeros(size(XALg,1),100);
for k=1:size(XALg,1)
    for z=1:100
    Xbucket(k,z)=sum(XALg(k,bucket_domain(z):bucket_domain(z+1)));
    end
end
    
for ind1=1:size(XALg,1)
    align{ind1}=ind1;
end

% Cluster the binned data 
    m=linkage(Xbucket,'weighted',distance_metric);
    figure, h=dendrogram(m,0);
    sample_order=str2num(get(gca,'XTickLabel'));
    tic

SegLength_range=[round(0.02/(ppm(2)-ppm(1))):(round(0.1/(ppm(2)-ppm(1))-round(0.02/(ppm(2)-ppm(1)))))/(100-1):round(0.1/(ppm(2)-ppm(1)))]; % range from 0.03 to 0.2, increase as less similar
MaxShift_range=[round(0.02/(ppm(2)-ppm(1))):(round(0.06/(ppm(2)-ppm(1))-round(0.02/(ppm(2)-ppm(1)))))/(100-1):round(0.06/(ppm(2)-ppm(1)))]; % should go from 0.01 to 0.03, increase as less similar
slack_range=[round(0.0015/(ppm(2)-ppm(1))):(round(0.01/(ppm(2)-ppm(1))-round(0.0015/(ppm(2)-ppm(1)))))/(100-1):round(0.01/(ppm(2)-ppm(1)))]; %range from 0.0015 to 0.01 ppm (5-30 points, at 64k on 20ppm range) should decrease as spectra get less similar 

for aligndex=1:size(m,1)
    
    XALg(isnan(XALg))=0;
    
    num1=m(aligndex,1);
    num2=m(aligndex,2);
    
    if round(100*m(aligndex,3))>=1
        spectra_distance=round(100*m(aligndex,3));
    else
        spectra_distance=1;
    end
    
    reference=mean([XALg(align{num1},:);XALg(align{num2},:)],1);
    SegLength=round(SegLength_range(spectra_distance));
    MaxShift=round(MaxShift_range(spectra_distance)); % should go from 0.005 to 0.03, increase as less similar
    slack=round(slack_range(spectra_distance)); %unsure about this parameter; let's work with it later, we think it should decrease
    NumSegs=round(size(X,2)/SegLength);
    switch upper(alignment_method);
       case('CCOW')
            XALg([align{num1},align{num2}],:)=CCOW([XALg(align{num1},:);XALg(align{num2},:)],reference,'SegLength',SegLength,'maxPeakShift',MaxShift,'Slack',slack);
       case('PARCCOW')
           XALg([align{num1},align{num2}],:)=parCCOW([XALg(align{num1},:);XALg(align{num2},:)],reference,'SegLength',SegLength,'maxPeakShift',MaxShift,'Slack',slack);
        case('ICOSHIFT')
            XALg([align{num1},align{num2}],:)=icoshift(reference,[XALg(align{num1},:);XALg(align{num2},:)],NumSegs,MaxShift);
        case('RAFFT')
            XALg([align{num1},align{num2}],:)=RAFFT([XALg(align{num1},:);XALg(align{num2},:)], reference, MaxShift, 1);
        case('PAFFT')
            XALg([align{num1},align{num2}],:)=PAFFT([XALg(align{num1},:);XALg(align{num2},:)], reference, SegLength, MaxShift);
    end
    
    align{ind1+1}=[align{num1},align{num2}];
    ind1=ind1+1;
    disp(['Completed ',num2str(aligndex),' of ',num2str(size(m,1))]);
    toc
    
end
end
