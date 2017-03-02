function [STOCSY_output] = STOCSY_mj_wrapper(X,ppm,targets,threshold,intT,figureOption,str)

%% Wrapper for STOCSY_mj_2

% This is a wrapper function for the STOCSY_mj_2 function. There is also
% a wrapperScript_STOCSY_mj_2 script to wrap this one. They should be
% kept in the same directory.
%{
% X:            spectral matrix (full spectral matrix is preferred for maximum
%               variance
% ppm:          ppm matrix, from low to high ppm values
% targets:      ppm values to be STOCSY'd (also called 'drivers')
% threshold:    STOCSY correlation threshold magnitude (0-1) for determining interactions.
%               STOCSY_mj_2 handles positive/negative correlations
% intT:         threshold for determining minimum peak height for counting correlation points as interactions 
% figureOption: 'generateFigures' pops out, saves (.fig), and closes a
%               customized figure, showing 1) Raw STOCSY output 2) Thresholded STOCSY 
%               3) All spectra with axes linked for easy investigation, for each driver in the run. 
%               This takes a lot more time and memory. Entering anything else, like 'noFigures', or
%               'enerateFigures' will suppress this and run MUCH faster.
%               Suppression recommended for runs when determining
%               thresholds, or just getting interactions. The function can
%               always be run on smaller lists of more interesting peaks as
%               primary drivers, if STOCSY figures are desired.
%addpath(genpath('')); % Add selected folders and subfolders to path
%}
%% Variables:
    % Wrapper variables
    % STOCSY_mj_2 variables
   
%%
     % Initialize Output Files
         % Interactions
                interactions_outFile = ['STOCSY_clusters_' str '_threshold_' num2str(threshold) '.csv'];          
                csvwrite(interactions_outFile,[]);
         % Network
                net_outFile = ['STOCSY_clusters_' str '_threshold_' num2str(threshold) '_binaryInteractions.csv'];
                csvwrite(net_outFile,[]);
         % Let us know where the Output is found and record in outputs
            fprintf(['Output Files for ' str ' stored in ','"',cd(),'"\n\n\n']); % so you know where the files are being put.

     % Initialize return structure
            STOCSY_output.step = str;
            STOCSY_output.targets = targets;
            STOCSY_output.threshold = threshold;
            STOCSY_output.interactions = [];
            STOCSY_output.interactions_outFile = interactions_outFile;
            STOCSY_output.X = X;
            STOCSY_output.ppm = ppm;
            STOCSY_output.intT = intT;
            STOCSY_output.figureOption = figureOption;
%% Get the data, Run the function
% Do STOCSY, thresholded and raw, generate and save figure, and return the
% data for the raw and thresholded STOCSYs and the responder peaklists in
% the "STOCSY_clusters_threshold_0.##.csv" file in the
% STOCSY_mj_Outputs_timestamp directory (above).
% target = targets(8);
    for i = 1:length(targets) % loop through all the targets in the list and run the function
        target = targets(i); % get the target from the list
        % Run function
            [~,~,interactions]=STOCSY_mj_2(target,X,ppm,threshold,intT,figureOption);
        % Get the driver in the front
            interactions = [targets(i)     interactions(interactions ~= targets(i))   ];
            
        % Print the result to file (difficult to store as a variable because of different row sizes)
            dlmwrite (interactions_outFile, interactions, '-append','precision',15);

    end
%% Get the results into the structure

STOCSY_output.interactions = csvread(interactions_outFile);

%% Reformat and write for Cytoscape (binary interactions)

interactionsList2cytoscape(csvread(interactions_outFile),net_outFile);
STOCSY_output.BinaryInteractions = csvread(net_outFile);

%% Calculate the correlation coefficiencts within all the interaction pairs
% {
    % There's a big issue here where rounding ppm and rounding the
    % interactions yields different numbers.
    index = zeros(size(STOCSY_output.BinaryInteractions));
    for i = 1:length(STOCSY_output.BinaryInteractions)
        [~,index(i,1)] = min(abs(ppm-STOCSY_output.BinaryInteractions(i,1)));
        [~,index(i,2)] = min(abs(ppm-STOCSY_output.BinaryInteractions(i,2)));
        a = corrcoef(X(:, index(i,1) ),X(:, index(i,2)));
        STOCSY_output.BinaryInteractions(i,3) = a(2);           
    end
%{
    try 
        for i = 1:13%length(STOCSY_output.BinaryInteractions)       
            a = corrcoef(X(:, find(ppmround==BIround(i,1)) ),X(:, find(ppmround==BIround(i,2)) ));
            STOCSY_output.BinaryInteractions(i,3) = a(2);           
        end
    catch
        fprintf('\n\n\nThe corr() function is apparently acting up. Check the BinaryInteractions field to make sure it worked.\n\n\n')
    end
%}
    % This error was being thrown because reading interactions in from file
    % resulted in zeros getting in (zero-filling for uneven rows). Issue 
    % comes up when ppm has no zero values in it, thus no index for them. Fixed by
    % using BinaryInteractions field, which has no zero-filling. 
    % Write to net_outFile
        csvwrite(net_outFile,STOCSY_output.BinaryInteractions);
%}
end

