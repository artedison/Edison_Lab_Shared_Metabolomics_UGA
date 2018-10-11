function [Sample] = initiateSamples(samples)
    Sample(1).ridges(1).RowInds = [];
    Sample(1).ridges(1).ColumnInds = [];
    Sample(1).ridges(1).WindowIndices = [];
    Sample(1).ridges(1).ppms = [];
    Sample(1).ridges(1).times = [];
    Sample(1).ridges(1).intensities = [];
    Sample(1).ridges(1).inputs = [];
    Sample(1).ridges(1).Adjusted = [];
    Sample(1).ridges(1).FilledGapInds = [];
    Sample(1).ridges(1).EndPointsAdded = [];
    
    if length(samples) > 1
        for i = 2:length(samples)
            Sample(i).ridges(1) =  Sample(1).ridges(1);
        end
    end
end