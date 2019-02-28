%% STEP 2: Ridge Tracing

% This is the workflow for manually guided ridge tracing of spectral data
% processed in STEP_1_processing_combine_samples.m . The algorithm works
% to cluster picked peaks in all spectra using their similarity as
% determined by their Euclidean coordinates in chemical shift (ppm), time,
% and intensity space. The user is then presented with a 3D figure
% showing all ridges (clusters of peaks). The user may click on a
% ridge to record it. Then, the currently selected ridges are adjusted.

% To run the script, the user:
% 1. Before running anything, ensure that the publicly available Edison Lab
% Metabolomics Toolbox and the data directories are added to the Matlab
% path. Each of the following sections are run as a block (e.g. CMND-Enter):
% 2. Run the "Setup" section once, and initiates the "Sample" storage
% structure. This will contain all relevant ridge information, including
% all parameters from the various steps in this workflow so that results
% can be reproduced. This is run at the beginning.

% 3. run the blocks "Smoothing and Peak Picking" through
% "Store current batch of ridges" for a given ROI in ALL samples. We
% were careful to pick ridges from different samples in the same order,
% and left empty placeholders where ridges in one sample did not exist in another.

% The "Smoothing and Peak Picking" step is run to optimize smoothing and
% peak-picking parameters. First, a region of interest (ROI) is defined
% in ppm units as a 2-element vector. The idea is that you pick ridges for
% a given region of interest, then pick the same ridges for all other
% samples. Next, a 2D Gaussian smoothing filter is applied:
% 'imgaussfilt(window,[1,1])', where [1,1] refers to the sigma for the
% Gaussian kernal in the [time,ppm] dimensions. The peakPickThreshold
% parameter can be any decimal value, positive or negative.
%
% The "Clustering" section operates on the smoothed matrix obtained in the
% preceding section for the ROI specified there. It involves manual tuning
% of clustering parameters, including:
% - numberOfRidges    % number of expected ridges. this needs to be relatively high for noisy regions
% - timeWeight        % weight in time domain. increasing this makes clusters reach across time
% - ppmWeight         % weight in ppm domain. increasing this makes clusters reach across the ppm dimension
% - intensityWeight   % weight in intensity domain. increasing this flattens the curves down to 2D (like the peakpick image)
% It is important to keep in mind that the latter three are all relative
% to each other and effectively define ratios, so their absolute values
% are arbitrary. The clustering results are displayed, and the user
% is prompted with a menu.
% Select "Pick Final Ridges" to record ridges by
% clicking on them. This works best by maximizing the figure window
% without rotation of the axes. When a ridge is selected, it will appear
% in bold. The user can also click a selected ridge to deselect it. When
% all ridges have been selected, press Enter. To close the figure without
% saving ridges, click "Cancel".
%
% The "Apply Ridge Corrections" section applies the following adjustments
% to the raw ridges:
% - remove side chains on ridges
% - use linear interpolation on gaps in ridges from smoothed data
% - map ridges to the original (non-smoothed) matrix
% - ensure all fields are sorted by timepoint
% - extend the ridge ends to fill the time series
%     (stored as a copy for optional use)
% - store all data in a structure with the original AND adjusted ridges
%
% Results of the corrections are displayed. The user must validate that
% The user may choose to adjust parameters in the corrections, or to omit
% different parts by modifying the code. The two relevant parameters are
% related to the mapping from smoothed to unsmoothed data:
% - windowWidth       % a positive integer indicating the number of data points on either side of the peak in the ridge after projection of a point onto unsmoothed data
% - viewWidth         % adjusts the number of ppms for viewing the window on either side of the peak
% First, the ridge point positions are projected directly onto the
% unsmoothed data. Then, a window of windowWidth data points are considered
% upfield and downfield of each ridge point. The maximum within each
% window is then taken as the new mapped ridge point. If the adjusted
% ridges are not satisfactory, the parameters should be adjusted.
%
% If they are, then the "Store current batch of ridges" section is used to simply
% concatenate the final ridges with the previously stored ones.
% Finally, the data must be saved in .mat format and passed to the next
% script (.m).
%
%   MJ OCT2018
%   YW 10/04/2018


%% 1. Setup
% load('sampleData.mat')
samples=[1,2,3,7,8,9]; % list of sample numbers from study
Samples=initiateSamples(samples);

%% 2. Smoothing and Peak Picking
sample=1; % index of "samples"
smoothingFilter='imgaussfilt(window,[1,10])';
peakPickThreshold=0;  % can be negative or positive; thresh = mean+(peakPickThreshold*sd)

% Run for different regions (get what you can out of them)
ppm=sampleData(samples(sample)).ppmR_1h1d;
matrix=sampleData(samples(sample)).Xcollapsed_1h1d;
times=sampleData(samples(sample)).timesCollapsed_1h1d;
%ROI = [5.45,5.7]; % in ppm units
ROI=[2.5,2.8];
%ROI = [2.0,2.5];

[output]=ridgeTracing_PeakPick1D(matrix,ppm,times,ROI,smoothingFilter,peakPickThreshold);

%% 3. Clustering
numberOfRidges=13; % this needs to be relatively high for noisy regions
timeWeight=1;     % increasing this makes clusters reach across time
ppmWeight=2;     % increasing this makes clusters reach across the ppm dimension
intensityWeight=500;  % increasing this flattens the curves down to 2D (like the peakpick image)
[newRidges]=ridgeTracing_clusterPeaks_interactive_2(output,numberOfRidges,timeWeight,ppmWeight,intensityWeight);

%% 4. Apply Ridge Corrections
matrix=sampleData(samples(sample)).Xcollapsed_1h1d;
currentppm=sampleData(samples(sample)).ppmR_1h1d;
windowWidth=5;
viewWidth=0.1; % in ppms around ridge extrema
[newRidges adjustedRidges]=ridgeCorrection(matrix,currentppm,newRidges,windowWidth,viewWidth,'plotFigs'); % 'plotFigs' to view them

%% Store current batch of ridges: comment one of the following linew
%% for the first one
Samples(sample).ridges=[newRidges];
%% for the afterwards ones
Samples(sample).ridges=[Samples(sample).ridges,newRidges];

%% Save the data
%save('ridges.mat')

%% Tracing Ridges
% See: STEP_3_combining_ridges.m
