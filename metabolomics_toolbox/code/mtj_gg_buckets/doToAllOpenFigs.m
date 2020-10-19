function doToAllOpenFigs(cmdStr)
%% doToAllOpenFigs

% Author: MTJ
% Version: 0.1
% Date: 2020
%
% Description:
%
%   Run the command (provided as a string) on all visible (open, not minimized or hidden) figures.
%   NOTE: when the command contains special characters (such as ' ), they
%   must follow MATLAB string interpretation standards. 
%
% Inputs:
%
%       cmdStr:  string to be evaluated for each figure 
%            (e.g.'set(gca,''xlim'',[1,2])'
%            NOTE: name-value params (e.g. 'xlim' above) will need extra
%            single quotes to be strings.
%
% Output:
%       
%       None; figures are modified according to the command (e.g. setting
%       xlims)
%
% Usage: 
%         
%         doToAllOpenFigs('set(gca,''xlim'',[7.5,9.2])')
%                 
% MTJ 2020
%% 



    figures = findall(0,'type','figure');
        for i = 1:length(figures)
           figure(i)
           eval(cmdStr);
        end    
end