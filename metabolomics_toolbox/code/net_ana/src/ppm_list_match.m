function [ppmmatch_ind ppmvec]=ppm_list_match(ppmlist,namelist,pattern,deltapm_threshold,except)
% The function will match the ppm in each cell of ppmlist and iteratively find the nearest ppm groups of ppm among all groups by pairwise distance. Only ppm has matched among all samples will remain. The matching will be carried out in filtered features (by pattern). The algorithm will iteratively compare the current 'reference' ppm with mean ppm from each replicate and include the matched ppm and update the 'reference' ppm adaptively.
%
% Arguments:
%         ppmlist: cell. numeric. the ppm list. Each cell is for one replicates and is a array (ntime*nfeature). must be provided
%         namelist: cell. string. the name list for each ppm. Each cell is for one replicates and contains names for all features. must be provided.
%         pattern: string. the pattern to search on namelist. default '[uU]nknown.*'
%         deltapm_threshold: float. ppm distance among estimated ppm and new ppm in the list. Default 0.01.
%         except: bool. whether exclude the pattern. default use the pattern instead of exclude it.
% Return:
%         ppmmatch_ind: the matched array for ppm in each cell (each column different ppm)
%         ppmvec: the estimated ppm value (by moving average among samples of different samples)
%
% Examples:
%
% ppmlist={[-0.000117598057077339	0.0573866054902868	0.876670179188786;
%           -0.000420251759958190	0.0573866054902868	0.876670179188786;
%           -0.000117598057077339	0.0573866054902868	0.876367525485905;
%           -0.000117598057077339	0.0573866054902868	0.876367525485905;
%           -0.000117598057077339	0.0573866054902868	0.876367525485905],
%          [-7.97819237602915e-05	-0.00552451331404291	0.0579973529059208;
%           -7.97819237602915e-05	-0.00552451331404291	0.0579973529059208;
%           -7.97819237602915e-05	-0.00552451331404291	0.0579973529059208;
%           -7.97819237602915e-05	-0.00552451331404291	0.0579973529059208;
%           -7.97819237602915e-05	-0.00552451331404291	0.0579973529059208],
%          [8.16588004817209e-05	-0.00536492014032353	0.0572707376789371;
%           -0.000220928918451879	-0.00536492014032353	0.0572707376789371;
%           -0.000220928918451879	-0.00536492014032353	0.0572707376789371;
%           -0.000220928918451879	-0.00536492014032353	0.0572707376789371;
%           8.16588004817209e-05	-0.00536492014032353	0.0572707376789371]};
% namelist={{'somethingelse'	'unknown'	'unknown'},{'somethingelse'	'unknown'	'unknown'},{'somethingelse'	'unknown'	'unknown'}};
% pattern='[uU]nknown.*';
% deltapm_threshold=0.01;
% [ppmmatch_ind ppmvec]=ppm_list_match(ppmlist,namelist,pattern,deltapm_threshold);
%
% Test:
% results = runtests('ppm_list_matchTest.m')
%
% Yue Wu 04/28/2020
% Tested with MATLAB R2018b

if ~exist('ppmlist','var')
  error('please provide input ppmlist cell');
end
if ~exist('namelist','var')
  error('please provide input namelist cell');
end
if ~exist('pattern','var')
  pattern='[uU]nknown.*';
end
if ~exist('deltapm_threshold','var')
  deltapm_threshold=0.01;
end
if ~exist('except','var')
  except=false;
end

ppmmatch_ind=[];
ppmmatch_ppm=[];
ppmvec=[];
nsample=length(namelist);
for isample=1:nsample
  namevec=namelist{isample};
  ppmmat=ppmlist{isample};
  matchind=cellfun(@(x) any(regexp(x, pattern)),namevec,'UniformOutput',true);
  if except
    matchind= ~matchind;
  end
  unanno_ind=find(matchind);
  ppmmean=mean(ppmmat,1);
  ppmmean_unannot=ppmmean(unanno_ind);%ppm mean vector for new replicate
  if isample>1
    % this block should be the same (logically) as ppm_vec_match().
    distmat=abs(pdist2(ppmvec',ppmmean_unannot'));
    %%%selecting out ppm points that are pairwise closest to each other
    [ppm_match_val1,ppm_match_ind1]=min(distmat,[],2);
    [ppm_match_val2,ppm_match_ind2]=min(distmat,[],1);
    ppm_match_val=[];
    ppm_match_ind=[];
    ppm_remain_ind=[];
    for ppm_min_i=1:length(ppm_match_ind1)
      if ppm_match_ind2(ppm_match_ind1(ppm_min_i))==ppm_min_i%check for pairwise match
        ppm_match_val=[ppm_match_val ppm_match_val1(ppm_min_i)];
        ppm_match_ind=[ppm_match_ind ppm_match_ind1(ppm_min_i)];
        ppm_remain_ind=[ppm_remain_ind ppm_min_i];
      end
    end
    thres_ind=find(ppm_match_val<deltapm_threshold);
    %%the two index to update record array and current ppm vector
    remainind=ppm_remain_ind(thres_ind);
    matchind=ppm_match_ind(thres_ind);
    ppmmatch_ind=[ppmmatch_ind(:,remainind); unanno_ind(matchind)];%%update the ppm index matrix
    ppmvec=ppmvec(remainind)*(isample-1)/isample+ppmmean_unannot(matchind)/isample;%%moving average mean ppm for estimating real ppm of the peak
    ppmmatch_ppm=[ppmmatch_ppm(:,remainind); ppmmean_unannot(matchind)];
  else
    ppmvec=ppmmean_unannot;
    ppmmatch_ind=unanno_ind;
    ppmmatch_ppm=ppmmean_unannot;
  end
  % whos ppmvec ppmmatch_ind ppmmean_unannot
end
