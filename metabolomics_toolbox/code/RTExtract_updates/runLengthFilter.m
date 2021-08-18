function [filteredVect,runbounds,runlengths,runinds] = runLengthFilter(vect,cutoff)

        % Get the boundaries, indices, and lengths of runs of true values
        % in a logical vector. 
        %
        % MTJ and YW 2020
        
%         cutoff = 3;
%         vect = [0 1 1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0];
        
        diffVect = diff([0; vect(:); 0]'); % Pad with zeros and take diff (with diff shift, this doesn't effect indices)
        runbounds = [find(diffVect == 1);find(diffVect == -1)-1]; % 1 means run starts on the current element (with diff shift), -1 means run ends on the preceding element
        runlengths = diff(runbounds,1)+1;
        runbounds(:,runlengths<cutoff) = []; % filter runs < cutoff length
        runlengths(runlengths<cutoff) = [];
        runinds = zeros(size(vect));
        for j = 1:size(runbounds,2)
            runinds(runbounds(1,j):runbounds(2,j)) = j;
        end
        filteredVect = ~(runinds==0);
end
