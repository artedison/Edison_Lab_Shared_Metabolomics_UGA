function [ant] = newAnt(matrix,noiseLevel,startPoint,halfWindow,varargin)

% Starting point

        % Default is random start point
%         if ~exist('startPoint','var') 
%             startPoint = round(rand(1) * numel(matrix)); % linear index of startpoint in matrix
%         end

        
        [m,n] = size(matrix);
        ant.r = zeros(m,1);
        ant.c = zeros(m,1);

        [r,c] = ind2sub([m,n], startPoint);
        
        ant.firstRow = r;       

    % Run loop for ant:
        
         ant.r(r) = r;
         ant.c(r) = c;
         r = r+1;
  
        if isempty(halfWindow)
            % Run deterministic
            % .01 s / 200 iterations
%             tic
            for i = r:m
                [ant.r(i),...
                 ant.c(i)] = ant_nextPos(ant.r(i-1),...     % changes, updates inside this fn
                                            ant.c(i-1),...     % changes
                                            n,matrix,...
                                            noiseLevel,...
                                            []); % half window for new point selection (empty = null)
            end
%             toc
            
        else
            % Run Stochastic
            rlist = rand(m,1); % pre-compute rand vals. outside of loop
            
%         tic
            for i = r:m
                [ant.r(i),...
                 ant.c(i)] = ant_nextPos(ant.r(i-1),...     % changes, updates inside this fn
                                            ant.c(i-1),...     % changes
                                            n,matrix,noiseLevel,... % stays the same each iteration
                                            halfWindow,rlist(i)); % half window for new point selection (stochastic)
            end
%         toc
        end
        inds = ant.r>0;
        ant.r = ant.r(inds);
        ant.c = ant.c(inds);
        
        ant.linearInds = sub2ind([m,n],ant.r,ant.c);
        ant.vals.inverted = matrix(ant.linearInds);
        ant.vals.reinverted = 1-ant.vals.inverted;
        
%         figure,hold on,
% %             plot(c,r)
%             set(gca,'xlim',[min(ppm),max(ppm)])
%             set(gca,'ylim',[1,m])
%             scatter(ppm(c),r,10,'r')        
%         scatter(ppm(ant(1).c(stpt:end)'),...
%                     ant(1).r(stpt:end),...
%                     10,'r')        
end
