function [optOB_out] = optimize_optBucket(X,ppm,bucketSizes,slacknessLevels)
%% optimize_optBucket
%
    % Author: Michael T. Judge
    % Version: 0.2
    % Tested on Matlab Version R2020b
    % Date: JUL2020
    %
    % Description:
    %       Applies opt_bucket without generating figures, with different 
    %       sb and sl parameter value combinations as provided, in order to 
    %       identify the optimal sb-sl combination for a given spectral 
    %       matrix. Each bucketSize and slacknessLevel combination is 
    %       tried. Then, filterBuckets_Peaks_opt() is used to filter out
    %       buckets without peaks. Results are then visualized with 
    %       plotOptBucket_optResult(), and the optimal parameter 
    %       combination is chosen using within the 
    %       plotOptBucket_optResult() function or set manually.
    %
    % Input:
    %       X                   spectral matrix
    %       ppm                 ppm vector matching matrix
    %       bucketSizes         vector of initial bucket sizes corresponding
    %                           to 'size_bucket' in opt_bucket()
    %       slacknessLevels     vector of slackness levels sizes corresponding
    %                           to 'slackness' in opt_bucket()
    %
    % Output:
    %       optOB_out           struct array containing results from each
    %                           opt_bucket run. Dimension 
    %                           length(bucketSizes)*length(slacknessLevels)
    %
    % Log:
    %
    % Example run:
    %
    %         size_bucket = 0.002:0.002:0.01;
    %         slackness = 0.45:0.13:0.99;
    % 
    %         [buckets] = optimize_optBucket(X,ppm,size_bucket,slackness);
    %         [buckets] = filterBuckets_Peaks_opt(ppm,buckets,peaks);
    % 

%% Generate buckets using a range of both params

%     sb = 0.002:0.002:0.01;
%     sl = 0.45:0.13:0.99;
    p = reportParams();
    ob_out = struct();
    
    for j = 1:length(bucketSizes)
            for i = 1:length(slacknessLevels) 
                n =  i + (j-1) * length(slacknessLevels);
                [ob_out(n).ZNN,ob_out(n).Z,ob_out(n).I_b,ob_out(n).S_b]=opt_bucket_noFigs(ppm,X,bucketSizes(j),slacknessLevels(i));
                    ob_out(n).sb = bucketSizes(j);
                    ob_out(n).sln = slacknessLevels(i);
            end
    end

    optOB_out.optimization.optParams_OB = p;
    optOB_out.optimization.results = ob_out;
    
end