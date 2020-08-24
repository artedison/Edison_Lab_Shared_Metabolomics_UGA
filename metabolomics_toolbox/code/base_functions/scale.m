function XS=scale(X,method,offset)

% XS=scale(X,method,offset)
%
% Normalize spectral matrix X to total area, single feature, or integral of
% set of features.
%
% Arguments:
%
% X               N x P matrix of spectral features for each sample 
% method          'log' for log2 fold-change vs. median scaling, 'logoff' 
%                 for offset log2 fold-change vs. median,
%                 'mc' for mean-centering (only), 'auto' for autoscaling
%                 (mean-center and univariance scaling), 'pareto' for
%                 pareto scaling. 'integral' for normalization to sum of
%                 set of features
% offset          optional feature for 'log' or 'logoff'.  'log' uses
%                 default offset of 0.0001 to prevent InF, 'logoff' uses automatically
%                 calculated median.
%
% Outputs: 
% XS              N x P matrix of scaled data
%
% MTJ edit 20200723: nans result from zero vals in several methods; replaced with
% zeros in result for robustness

switch method
    case 'mc'
        XS=X-repmat(mean(X),[size(X,1),1]);
        
    case 'auto'
        XMC=X-repmat(mean(X),[size(X,1),1]);
        XS=XMC./repmat(std(X),[size(X,1),1]);
    case 'range'
        XMC=X-repmat(mean(X),[size(X,1),1]);
        Xmin=repmat(min(X,[],1),[size(X,1),1]);
        Xmax=repmat(max(X,[],1),[size(X,1),1]);
        XS=XMC./(Xmax-Xmin);
    case 'pareto'
        XMC=X-repmat(mean(X),[size(X,1),1]);
        XS=XMC./repmat(sqrt(std(X)),[size(X,1),1]);
    case 'vast'
        XMC=X-repmat(mean(X),[size(X,1),1]);
        XS=XMC./repmat((std(X)),[size(X,1),1]);
        XA=repmat(mean(X),[size(X,1),1])./repmat((std(X)),[size(X,1),1]);
        XS=XS.*XA;
    case 'level'
        XMC=X-repmat(mean(X),[size(X,1),1]);
        XS=XMC./repmat(mean(X),[size(X,1),1]);
        
    case 'log'
        X=abs(X);
        if exist('offset')~=1
            offset=0.0001;
        end
        median_X=median(X+offset);
        FC=zeros(size(X));
        for k=1:size(X,1)
            FC(k,:)=[X(k,:)+offset]./median_X;
        end
        XS=log2(FC);
    case 'logoff'
        X=abs(X);
        if exist('offset')~=1
            offset=median(X(X>0));
        end
        median_X=median(X+offset);
        FC=zeros(size(X));
        for k=1:size(X,1)
            FC(k,:)=[X(k,:)+offset]./median_X;
        end
        XS=log2(FC);
    case 'power'
        X=abs(X);
        if exist('offset')~=1
            offset=0.0001;
        end
        median_X=median(X+offset);
        FC=zeros(size(X));
        for k=1:size(X,1)
            FC(k,:)=[X(k,:)+offset]./median_X;
        end
        XS=pow2(FC);
    case 'poweroff'
         X=abs(X);
        if exist('offset')~=1
            offset=median(X(X>0));
        end
        median_X=median(X+offset);
        FC=zeros(size(X));
        for k=1:size(X,1)
            FC(k,:)=[X(k,:)+offset]./median_X;
        end
        XS=pow2(FC);
    
    case 'log10'
        X=abs(X);
       XMC=log(X);
       XS=XMC-(repmat(mean(X),[size(X,1),1]));
    case 'pwr'
        X=abs(X);
        XMC=sqrt(X);
        XS=XMC-(repmat(mean(X),[size(X,1),1]));
    case 'optim_glog'
        X=abs(X);
        [XS,~,~]=optimise_glog_parameter(X);
        XS=XS';
    case 'glog'
        X=abs(X);
        [XS,~,~]=glog_opt(X);

end

        XS(isnan(XS)) = 0; % MTJ edit 20200723: nans result from zero divisors, which happens frequently. 

end