%% STEP 3: Combining Ridges

% Notes:
    % It seems best to not calculate mean and SD across samples if
    % their timepoints aren't the same. Maybe we can do this in a 2D way,
    % but that conflates P0 and biological error.

    % It seems best to sum all ridges annotated as a given compound within
    % a sample, as these truly have the same timepoints and this will
    % reduce the effect of noise and baseline. With this, each compound
    % needs to be reported with the number of ridges used to quantify it.
    % The caveat here is that the ridges are all different lengths. As
    % such, summing is not really possible, and statistics such as the
    % number of ridges per compound are not entirely helpful.

    % The Adjusted.ExtendedIntensities field should not be used, as it is
    % complicated to access. Instead, these data (and more) are now also
    % stored in Samples.adjustedRidges.

% load('ridges.mat')
Samples=Samples2;
samples=[1,2,3,7,8,9]; % maps j to correct index in sampleData

%% 1. Adjust the time vectors according to inoculation-start time lag
% Sample #, time lag (days), sampleData index
lags = [4   0.013194444 1
        5   0.011805556 2
        6   0.011805556 3
        10  0.018055556 7
        11  0.01875 8
        12  0.015277778 9];
lags(:,2)=lags(:,2)*24; % convert days to hours

% For each sample correct time
for s=1:length(Samples)
    for r=find(~cellfun(@isempty,{Samples(s).adjustedRidges.LinearInds})) % only loop through actual ridges
        Samples(s).adjustedRidges(r).Times=Samples(s).adjustedRidges(r).Times+lags(s,2);
    end
end

%% 2. Set up the ridge-compound peaks map. Assign each ridge a name and compound.
sampleKey={'aerobic','aerobic','aerobic','anaerobic','anaerobic','anaerobic'};
ridgeNames={'unknown-179','valine 1' 'valine 2' 'valine 3' 'valine 4' 'isoleucine 1' 'isoleucine 2' 'alanine 1' 'alanine 2' 'lactate 1' 'lactate 2' 'arginine 1' 'arginine 2' 'arginine 3' 'arginine 4' 'arginine 5' 'arginine 6' 'arginine 7' 'arginine 8' 'uracil 1' 'uracil 2' 'glucose-1-phosphate 1' 'glucose-1-phosphate 2' 'glucose-1-phosphate 3' 'glucose-1-phosphate 4' 'adenosine 3' 'adenosine 4' 'tyrosine 1' 'tyrosine 2' 'phenylalanine 1' 'phenylalanine 2' 'uridine 1' 'uridine 2' 'choline' 'trehalose 1' 'trehalose 2' 'guanosine' 'citrate 1' 'citrate 2' 'citrate 3' 'citrate 4' 'succinate' 'glutamate' 'fumarate' 'formate' 'glucose 1' 'glucose 2'  'glucose 3' 'glucose 4'  'glucose 5' 'glucose 6' 'glucose 7' 'glucose 8' 'ethanol 1' 'ethanol 2' 'ethanol 3' 'ethanol 4' 'ethanol 5' 'ethanol 6' 'ethanol 7' };
ind=[179 75 76 79 80 77 78 91 92 86 87 95 97 98 99 100 101 102 103 55 56 64 65 66 67 14 16 41 42 32 33 23 24 151 152 153 154 158 159 160 161 162 163 164 13 166 167 168 169 170 171 156 157 172 173 174 175 176 177 178];
annotations=regexprep(ridgeNames,'( )\d',''); % get rid of the peak number for compounds with multiple peaks
compoundList=unique(annotations); % how many/which compounds have we quantified?

% numberOfProtons=repmat(1,1,length(annotations));  % this is a dummy vector for number of protons represented by each peak (this is actually not applied here)
unknowns=setdiff(1:length(Samples(1).ridges),ind);
fullInd=[ind,unknowns];
ridgeNames=[ridgeNames,num2cell(unknowns)];
% numberOfProtons=[numberOfProtons,repmat(1,1,length(unknowns))];
annotations=[annotations,repmat({'unknown'},1,length(unknowns))];

% To append a cell array, loop through Samples and ridges
[fullInd,sortInds]=sort(fullInd);
ridgeNames=ridgeNames(sortInds);
% numberOfProtons=numberOfProtons(sortInds);
annotations=annotations(sortInds);
for j=1:length(Samples)
    for i=1:length(Samples(j).ridges)
        Samples(j).adjustedRidges(i).RidgeName = ridgeNames{i};
        Samples(j).adjustedRidges(i).Annotation = annotations{i};
    end
end

%% Sum all ridges (normalized to #protons) for each compound (excluding unknowns)
[compounds]=combineRidges(Samples,compoundList,samples,sampleKey);

%% Add plot colors to the structure
% Make the colors
%uisetcolor() % Used this to get the desired colors
blue=[0,.3,1];
yellow=[1    0.7490    0.0510];%
red=[0.8,0,0];

% Add to the compounds structure
for c = 1:length(compounds)
    if strcmp(compounds(c).Condition,'aerobic')
        compounds(c).PlotColor = red;
    else
        compounds(c).PlotColor = blue;
    end
end

%% Save the data
%save('compounds.mat')

%% Plotting Compound dating
% See: STEP_4_plotting.m
