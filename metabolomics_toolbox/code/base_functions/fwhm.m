function [Fullw] = fwhm(X,ppm,shift1,shift2,label)

% Format XTitles -- Ignore if no input
if nargin == 5
    % change row vectors to column vectors
    if size(label,1) == 1
        label = label';
    end
else
    label = [];
end

ppm_range = ppm;

% Remove regions on ppm_range
for i = 1:length(ppm_range)
    if ppm_range(i) < shift1 || ppm_range(i) > shift2
        ppm_range(i) = 0;
    end
end

% find max of each spectra
for i = 1:size(X,1)
    data = X(i,:);      % current spectra
    data_range = data(logical(ppm_range));    % find only peak b/n shift1 and shift2
    peak = max(data_range);	% peak of spectra
    half_peak = peak/2; % half peak
    index_of_peak = find(data==peak);
    
    % index_of_half_peak consists of the 1st and 2nd ppm at half peak
    index_of_half_peak = [index_of_peak, index_of_peak];   % starting point
    
    % find 1st index_at_half_peak -- downstream of peak
    while half_peak < data(index_of_half_peak(1))
        index_of_half_peak(1) = index_of_half_peak(1) - 1;
    end
    
    % find 2nd index_at_half_peak -- upstream of peak
    while half_peak < data(index_of_half_peak(2))
        index_of_half_peak(2) = index_of_half_peak(2) + 1;
    end
    
    % ppm at half peaks
    ppm1 = ppm(index_of_half_peak(1));
    ppm2 = ppm(index_of_half_peak(2));
    
    % full width
    Fullw(i,1) = ppm2 - ppm1;  
end

% Add label if exists
if length(label)
    % convert to cell
    if iscell(label)
        temp = Fullw;
        Fullw = cell(length(temp),1);
        for i = 1:length(temp)
            Fullw{i} = temp(i);
        end
    end
    % concatenate label and full width
    Fullw = [Fullw,label];
end



    