function [finX,finT] = subset_X(X,T,sel,Yvec)

% Author: Rahil Taujale
% Version: 0.1
% Date: 01/30/2017

% Description:
%       Selects a subset of X based on a provided list of Yvec or sample
%       IDs or sample group names.
%
% Input: 
%       X       : stack of 1D spectra
%       T       : Table with the sample information
%                   5 columns with following exact headers:
%                    Run_ID,Sample_ID,Sample_desc,Sample_grp,Yvec
%       sel     : 3 options.
%           1) A cell array with a list of Sample_grp names to be selected,
%           corresponds to the Sample_grp column in table T.
%               Eg: sel={'L1','L2'}
%           2) A numeric vector with a list of Yvec to be selected
%               Eg: sel=[5 6 7]
%           3) A cell array with a list of Sample_ID to be selected.
%       Yvec    : 'Y' if a list of Yvec is provided
%                 'S' if a list of Sample_ID is provided
%                 '' leave unassigned if Sample_grp is provided
%
% Output: 
%       A new matrix with added data points.
%       A new table T to be ised with htis new data matrix X
%
% Example run: [newSubX,newSubT]=subset_X(JointX,T,selected_sets,'Yvec');

    [M,N]=size(sel);
    finX=[];
    finT=[];
    for f=1:N
        if strncmpi(Yvec, 'Yvec',1)
            Xsel=X(T.Run_ID(T.Yvec==sel(f)),:);
            Tsel=T(T.Run_ID(T.Yvec==sel(f)),:);
        elseif strncmpi(Yvec, 'sample',1)
            Xsel=X(T.Run_ID==sel(f),:);
            Tsel=T(T.Run_ID==sel(f),:);
        else
            Xsel=X(T.Run_ID(strcmp(T.Sample_grp, sel(f))),:);
            Tsel=T(T.Run_ID(strcmp(T.Sample_grp, sel(f))),:);
        end
        finX=vertcat(finX,Xsel);
        finT=vertcat(finT,Tsel);
    end
    finT.Run_IDOrig=finT.Run_ID;
    finT.Run_ID=(1:1:length(finT.Run_IDOrig))';
end
    