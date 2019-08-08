function [compounds] = combineRidges_new(Samples,compoundList,samples,sampleKey)
%% this function is adpated from combineRidges.m and changed to make it works for new data structure
%% Sample: the ridge structure of the ridge tracing method
%% compoundList: list of name of compound to be quantified
%% samples: the sample index
%% sampleKey: sample name
%% michael judge
%% yue wu

% Pre-allocate the struct
numSamples=length(Samples);
numCmpds=length(compoundList);
numCombs=numCmpds*numSamples;
compoundList=sort(compoundList);
compounds=struct('Name',sort(repmat(compoundList',numSamples,1)),...
        'SampleNumber',num2cell(repmat(samples',numCmpds,1)),...
        'Condition',repmat(sampleKey',numCmpds,1),...
        'PlotColor',repmat(sampleKey',numCmpds,1),...
        'Times',cell(numCombs,1),...
        'RidgeIndices',cell(numCombs,1),...
        'RowIndices',cell(numCombs,1),...
        'ColumnIndices',cell(numCombs,1),...
        'ppms',cell(numCombs,1),...
        'Intensities',cell(numCombs,1),...
        'trajectory_RidgeIntensities',cell(numCombs,1),...
        'trajectory_ScaledRidgeIntensities',cell(numCombs,1),...
        'AverageTrajectory_times',cell(numCombs,1),...
        'AverageTrajectory_std_raw',cell(numCombs,1),...
        'AverageTrajectory_intensities_raw',cell(numCombs,1),...
        'AverageTrajectory_intensities_shifted',cell(numCombs,1),...
        'AverageTrajectory_intensities_scaled',cell(numCombs,1),...
        'AverageTrajectory_overlapRows',cell(numCombs,1),...
        'AverageTrajectory_ridgesUsed',cell(numCombs,1),...
        'AverageTrajectory_numberOfRidges',cell(numCombs,1));

for c=1:length(compoundList)
  for s=1:length(Samples)
      % Row indices for this compound-sample combination
      ind=s+numSamples*(c-1);
      % For each sample, we need to search the list of all compounds
      datasamp=Samples(s).ridges;
      compnamdvec=[];
      datasamp(1)=[];
      for i=1:length(datasamp)
        compnamdvec=[compnamdvec unique(datasamp(i).result.names)];
      end
      compoundInds=find(strcmp(compnamdvec,compoundList{c}));
      % Only go through this process if there is actually a ridge:
      if ~isempty(compoundInds)
          % transfer the data to the compounds structure:
          % CHECK TO ENSURE THAT FIELDS end up the same size!!! (e.g. 1
          % column of cells)
          for compoundInds_i=1:length(compoundInds)
            compoundInds_ele=compoundInds(compoundInds_i);
            compounds(ind).Times{compoundInds_i}=datasamp(compoundInds_ele).result.time';
            % compounds(ind).RidgeIndices=compoundInds;
            compounds(ind).RowIndices{compoundInds_i}=datasamp(compoundInds_ele).result.rowind';
            compounds(ind).ColumnIndices{compoundInds_i}=datasamp(compoundInds_ele).result.colind';
            compounds(ind).ppms{compoundInds_i}=datasamp(compoundInds_ele).result.ppm';
            intensityhere=datasamp(compoundInds_ele).result.intensity';
            minintensityhere=min(intensityhere);
            if minintensityhere<0
              intensityhere=intensityhere-minintensityhere;
            end
            compounds(ind).Intensities{compoundInds_i}=intensityhere;
          end
          compounds(ind).Times=compounds(ind).Times';
          compounds(ind).RowIndices=compounds(ind).RowIndices';
          compounds(ind).ColumnIndices=compounds(ind).ColumnIndices';
          compounds(ind).ppms=compounds(ind).ppms';
          compounds(ind).Intensities=compounds(ind).Intensities';
          timeRange=unique([compounds(ind).Times{:}]); % these are the times for which there are one or more ridge points. unique() sorts.
          compounds(ind).AverageTrajectory_times=timeRange;
          % Do the averaging
          compounds(ind).AverageTrajectory_intensities_shifted=zeros(1,numel(timeRange));
          compounds(ind).AverageTrajectory_intensities_scaled=zeros(1,numel(timeRange));
          compounds(ind).AverageTrajectory_intensities_raw=zeros(1,numel(timeRange));
          compounds(ind).AverageTrajectory_std=zeros(1,numel(timeRange));
          compounds(ind).AverageTrajectory_numberOfRidges=zeros(1,numel(timeRange));
          compounds(ind).trajectory_RidgeIntensities=cell(1,numel(timeRange));
          compounds(ind).AverageTrajectory_ridgesUsed=cell(1,numel(timeRange));
          % Set up to scale each ridge by the difference between it and the max for those indices
          % Starting with the largest ridge, calculate the linear
          % indices of the regions of overlap with other ridges.
          % Assuming that the points of highest overlap are the
          % highest points of the ridges (which should track each
          % other), this region should be used to make the adjustment
          % (scaling). This is because this point will be where the
          % smallest peaks actually get strong enough to see.
          % Find the region
          rowInds=[compounds(ind).RowIndices{:}];
          counts=sum(rowInds==rowInds');
          overlapInds=find(counts==max(counts)); % this is the region used to calculate the means; the region of highest overlap
          overlapRows=unique(rowInds(overlapInds)); % these are the timepoints with the most overlap across ridges
          % Get the inds of that region in each cell of RowIndices
          regionIndsByCell=cell(length(compoundInds),1);
          for r=1:length(compoundInds)
              [~,regionIndsByCell{r}]=intersect(compounds(ind).RowIndices{r},overlapRows);
          end
          % Get the means of that region for each ridge
          means=zeros(length(compoundInds),1);
          for r=1:length(compoundInds)
              means(r)=mean(compounds(ind).Intensities{r}(regionIndsByCell{r}));
              %ranges(r) = range(compounds(ind).Intensities{r}(regionIndsByCell{r}));
          end
          % Apply the shift (shift the ridges by addition)
          magnitudeShifts=max(means)-means;
          for r=1:length(compoundInds)
              compounds(ind).trajectory_ScaledRidgeIntensities_shifted{r}=compounds(ind).Intensities{r}+magnitudeShifts(r);
          end
          % Apply the shift (scale the ridges by ratio of means)
          magnitudeShifts=max(means)./means;
          for r=1:length(compoundInds)
              compounds(ind).trajectory_ScaledRidgeIntensities_scaled{r}=compounds(ind).Intensities{r}*magnitudeShifts(r);
          end
          % Calculate the mean of the scaled ridges at each timepoint
          % Get the intensities in a linear array
          intensities=cell2mat(compounds(ind).Intensities'); % THIS REQUIRES THAT THE FIELDS ARE THE SAME SIZES
          intensitiesShifted=cell2mat(compounds(ind).trajectory_ScaledRidgeIntensities_shifted); % THIS REQUIRES THAT THE FIELDS ARE THE SAME SIZES
          intensitiesScaled=cell2mat(compounds(ind).trajectory_ScaledRidgeIntensities_scaled); % THIS REQUIRES THAT THE FIELDS ARE THE SAME SIZES
          % Loop through the timepoints and calculate the means
          for t=1:numel(timeRange)
              % Find all points that belong to the current time. We'll do this with unlisted cell arrays
              % (given that they have equivalent sizes) and linear indexing.
              thisTimeInds=find(cell2mat(compounds(ind).Times')==timeRange(t)); % unlists them cell by cell
              compounds(ind).trajectory_RidgeIntensities{t}=intensities(thisTimeInds);
              % Raw intensities average and sd()
              compounds(ind).AverageTrajectory_intensities_raw(t)=mean(intensities(thisTimeInds));
              compounds(ind).AverageTrajectory_std_raw(t)=std(intensities(thisTimeInds));
              % Shifted intensities average and sd()
              compounds(ind).AverageTrajectory_intensities_shifted(t)=mean(intensitiesShifted(thisTimeInds));
              compounds(ind).AverageTrajectory_std_shifted(t)=std(intensitiesShifted(thisTimeInds));
              % Scaled intensities average and sd()
              compounds(ind).AverageTrajectory_intensities_scaled(t)=mean(intensitiesScaled(thisTimeInds));
              compounds(ind).AverageTrajectory_std_scaled(t)=std(intensitiesScaled(thisTimeInds));
              % Other info:
              compounds(ind).AverageTrajectory_numberOfRidges(t)=numel(thisTimeInds);
              compounds(ind).AverageTrajectory_ridgesUsed{t}=compoundInds(find(cell2mat(cellfun(@any,cellfun(@(x) x==timeRange(t), compounds(ind).Times, 'UniformOutput', 0),'UniformOutput', 0))));
              compounds(ind).AverageTrajectory_overlapRows=overlapRows;
          end
      end
  end
end
compounds(find(cellfun(@isempty,{compounds.Times})))=[];  % delete empty rows
end
