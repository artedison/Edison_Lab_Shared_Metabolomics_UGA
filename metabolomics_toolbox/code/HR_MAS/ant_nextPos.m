function [newr,newc] = ant_nextPos(r,c,n,matrix,noiseLevel,halfWindow,rnum)

%% One step       
    % Can move down a row in the matrix, or left and right
        % Define the region locally:
        %   na  na na
        %     - - - 
        %   l * p * r
        %     * * *
        %   ll  b  lr
        
            % matrix lininds of l, ll, b, r, lr:
 
%                 locreg = [p-m;...
%                           p-m+1;...
%                           p+1;...
%                           p+m;...
%                           p-m+1];

%% Shake up start point:

    % If halfWindow is provided, add random jump behavior.
    % For vis/troubleshooting, commented plot commands available. If
    % halfWindow is used, rnum must be provided as well (rand(1),
    % pre-computed outside of the loop).
   
        if ~isempty(halfWindow)
                % Shake things up: read comments for details. Basically, calculate a
                % random alternate start column jumpPt within window around current 
                % start column c. Walking from c to jumpPt, assess if any local extrema 
                % need to be traversed where a valley -> peak transition > noiseLevel.
                % If not, accept the new point and set c = jumpPt.
            %
%                                         figure,hold on,
%                                         plot(matrix(r,:)') %,r+1
%                                         scatter(c,matrix(r,c),'r*')

                        % Define acceptable window in which to look for jump points
                         window = max([1,c-halfWindow]) : min([c+halfWindow,size(matrix,2)]);
             
%                                         xline([window(1),max(window)])
%                                         xline(c,'r')

                            % Randomly choose a jump point within window
%                                 jumpPt = window(max([1,round(rand(1) * length(window))]));
                                jumpPt = window(max([1,round(rnum * length(window))]));

%                                         plot(jumpPt,matrix(r,jumpPt),'*r')
%                                         xline(jumpPt,'--r')

                            % Determine directionality
                                if jumpPt>c
                                    reg = [c,jumpPt];
                                    flip = 1;
                                else
                                    reg = [jumpPt,c];
                                    flip = -1; % this will effectively reverse the direction of the diff later
                                end

                            % Get region between current point c and proposed jump point    
                                regvals = matrix(r, reg(1):reg(2));

                                if length(regvals)>3
                                
                                    [pks,plocs] = findpeaks(regvals);
                                        plocs = reg(1) + plocs - 1; %+ flip?
            %
%                                             plot(plocs,pks,'*k')

                                    [vals,vlocs] = findpeaks(1-regvals); % valleys to peaks
                                        vals = 1-vals; % peaks to valleys
                                        vlocs = reg(1) + vlocs - 1;
            %                            
%                                         plot(vlocs,vals,'*b')
                                else
                                    plocs = [];
                                    vlocs = [];
                                    pks = [];
                                    vals = [];
                                end


                            % Including c and jumpPt, determine if any valley -> peak
                            %   height diff is > noiseLevel. If so, don't jump.
                                [~,inds] = sort([reg(1),plocs,vlocs,reg(2)]);
                                extrema = [regvals(1),pks,vals,regvals(end)];

                                % Do the comparison. flip = -1 reverses order of
                                % diff() from l2r -> r2l
                                if ~any(flip*diff(extrema(inds)) > noiseLevel) % moving down OK, not up

                                    c = jumpPt; % make the jump; set c (start column) to jumpPt
            %                        
                                        % If jump is made, circle the new point.
%                                         plot(c,matrix(r,c),'or','MarkerSize',15)
                                end
                                
        end
        
        
%% Find local minimum on current row:

            rowReg = [c-1,c,c+1];
                    rowReg = rowReg(rowReg <= n & rowReg >= 1); % make sure it doesn't walk off the side
                [~,ind] = min(matrix(r,rowReg));
                localMin = rowReg(ind);
                
            while localMin ~= c
                c = localMin;
                rowReg = [c-1,c,c+1];
                    rowReg = rowReg(rowReg <= n & rowReg >= 1); % make sure it doesn't walk off the side
                [~,ind] = min(matrix(r,rowReg));
                localMin = rowReg(ind);
            end
     
%             scatter(newc,r,10,'r')
            
%% Move forward left, forward straight, or forward right
% By subscripts:
                            
        % Get vals, add noise, calc movements are calculated by assessing 
        % the landscape with added noise (or smoothing):
        
%             rows = [r+1,r+1,r+1];
            newr = r+1;         % this is mandatory
            
            cols = [localMin - 1, localMin, localMin+1];    % look at local neighborhood in next row
            cols = cols(cols <= n & cols >= 1); % make sure it doesn't walk off the side
            
            [~,minreg] = min(matrix(newr,cols));
            newc = cols(minreg);

%             scatter(newc,r+1,10,'r')
% Linearized:
%             p = sub2ind([m,n],m,lowerpoint);
%                             
%         % Get vals, add noise, calc movements are calculated by assessing 
%         % the landscape with added noise (or smoothing):
%             locreg = [p-m+1,p+1,p+m+1];                             % take min of
% %             [~,minreg] = min(matrix(locreg) + ...                 % vals + 
% %                             rand(1,3)*noiseLevel*2 - noiseLevel); % noise
%             [~,minreg] = min(matrix(locreg)); % no noise
%                       
%             nextp = locreg(minreg); % get location
        
end
