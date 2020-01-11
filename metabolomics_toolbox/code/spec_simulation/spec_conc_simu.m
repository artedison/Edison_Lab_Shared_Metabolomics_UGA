function [simures] = spec_conc_simu(args)
%% the function to generate concentration matrix (sample X compound)
% arguments: args
%     args.conc_range: the whole dynamic range of all compounds. default [0 1000]
%     args.nsample: number of samples. default 50
%     args.simutype: type of simulation: 'random_chemical_sample' OR 'biological_sample'. One of the type must be provided.
%%%         'random_chemical_sample': each sample is considerred independent and compound concentrations are random distributed. For each compound, a random range was chosen and the range is used to transform sample from normal distribution to produce concentrations for each time point.
%%%         'biological_sample': N biological samples with some compound concentrations different from the samples
%     args.factor_std: the relative variance for compound concentration in samples. default 0.3
%     args.compd_vec: the simulated compound list. Must be provided. {'A' 'B' 'C'}
%     args.modif_conc_list: concentration that are changed between groups. For each compounds {compound, relative mean concentration in different samples}
%         EX:  cell2table({'L-arginine' 1 3; 'Choline' 2 1});
%     args.corr_conc_list: correlation between compound and phenotypes. {group index, compound, phenotype index, simulated correlation}
%         EX:  cell2table({2 'Glycine' 1 0.8});
%     args.sampleindvec: sample index vector. EX. [1 1 1 2 2 2] two groups and each with 3 samples
%     args.pheno_mat: the phenotype matrix (sample X phenotypes). Should be provided if corr_conc_list is not empty
%
% Return: simures:
%               conc_mat: the concentration matrix (sample X compound)
% YUE WU 12112019

if ~isfield(args,'conc_range')
  args.conc_range=[0 1000];
end
if ~isfield(args,'nsample')
  args.nsample=50;
end
if ~isfield(args,'simutype')
  error('please provide simulation types');
end
if ~isfield(args,'factor_std')
  args.factor_std=0.3;
end
if ~isfield(args,'compd_vec')
  error('please provide compound list');
end
if ~isfield(args,'modif_conc_list')
  args.modif_conc_list=cell2table({nan nan nan;});
end
if ~isfield(args,'pheno_mat') && ~isfield(args,'corr_conc_list')
  error('please provide phenotype matrix');
end
if ~isfield(args,'corr_conc_list')
  args.corr_conc_list=cell2table({nan nan nan nan;});
end
if ~isfield(args,'sampleindvec')
  error('please provide sample vector');
end

conc_range=args.conc_range;
nsample=args.nsample;
simutype=args.simutype;
factor_std=args.factor_std;
compd_vec=args.compd_vec;
modif_conc_list=args.modif_conc_list;
corr_conc_list=args.corr_conc_list;
sampleindvec=args.sampleindvec;
pheno_mat=args.pheno_mat;

factor_std_vec=[repmat(factor_std,1,10) repmat(factor_std*1.5,1,10)];
%%% simulating different compound concentration for different samples
conc_mat=zeros(length(sampleindvec),length(compd_vec));
for compdi=1:length(compd_vec)
  compd=compd_vec{compdi};
  %%random generate range for the compound between a preset range
  crang=sort(rand(1,2)*(conc_range(2)-conc_range(1))+conc_range(1));%%defined the range for a compound (default)
  if strcmp(simutype,'random_chemical_sample')
    for samplei=1:length(sampleindvec)
      %%randomness in existence * randome sampled concentration
      cexist=randi(2)-1;
      cconc=exp(normrnd(0,factor_std_vec(compdi)))*(crang(2)-crang(1))+crang(1);%%log normal distributed for a small factor_std it is close to linear
      conc_mat(samplei,compdi)=cexist*cconc;
    end
  elseif strcmp(simutype,'biological_sample')
    modif_compd_ind=find(strcmp(table2cell(modif_conc_list(:,1)),compd));
    corr_compd_ind=find(strcmp(table2cell(corr_conc_list(:,2)),compd));
    %%count samples
    sampleuniq=sort(unique(sampleindvec));
    sampcount=[];
    for sampleuniqi=1:length(sampleuniq)
      sampleuniqele=sampleuniq(sampleuniqi);
      sampcount=[sampcount length(find(sampleindvec==sampleuniqele))];
    end
    %%group differences
    ratio=[];
    if length(modif_compd_ind)~=0
      modif_compd_vec=table2cell(modif_conc_list(modif_compd_ind,:));
      ratio=[];
      for sampleuniqi=1:length(sampleuniq)
        ratio=[ratio repmat(modif_compd_vec{1+sampleuniqi},1,sampcount(sampleuniqi))];
      end
    else
      ratio=repmat(1,1,length(sampleindvec));
    end
    %%exp to make sure positive
    conc_mat(:,compdi)=(exp(normrnd(0,factor_std_vec(compdi),length(sampleindvec),1))*(crang(2)-crang(1))+crang(1)).*ratio';
    %%correlations between compound and phenotypes
    if length(corr_compd_ind)~=0
      corr_conc_vec=table2cell(corr_conc_list(corr_compd_ind,:));
      corrind=find(sampleindvec==corr_conc_vec{1});
      incorrind=find(sampleindvec~=corr_conc_vec{1});
      xvec=pheno_mat(corrind,corr_conc_vec{3});
      corrval=corr_conc_vec{4};
      ypre=conc_mat(corrind,compdi);
      %%%% adpated from stackexchange: CrossValidation: whuber 2019
      % xvec=rand(50,1);
      % ypre=rand(50,1);
      % corrval=0.8;
      md=fitlm(xvec,ypre);
      xg=md.Residuals.Raw;
      ysimu=corrval*std(xg)*xvec+xg*std(xvec)*sqrt(1-corrval^2);
      ysimu=ysimu-mean(ysimu)+mean(conc_mat(incorrind,compdi));
      % corr(xvec,ysimu);
      %%%%
      conc_mat(corrind,compdi)=ysimu;
    end
  end
end
simures=struct();
simures.conc_mat=conc_mat;
