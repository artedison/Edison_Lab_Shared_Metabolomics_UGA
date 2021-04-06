function [cinfo] = customColormap(x,varargin)

%% customColormap
%
%   Generates continuous rgb colors and custom colormap for a vector/matrix x
%   Linearly interpolated between provided or default colors for true color
%   correspondence. R, G, and B values are interpolated separately, then
%   re-combined.
% 
% Inputs: 
%     
%     x           vector of data values (e.g. measurements)
%     varargin    
%         'colors'        n x 3 matrix of RGB values defining n colors between which to
%                         build the map. Default: 
%                         Blue (low vals), White (middle vals), Red (high vals)
%         'setPoints'     values to base the color positions on (data values where the colormap
%                         is equal to the corresponding color in 'colors'. Default:
%                         [min(x),median(x), max(x)].
%                         * NOTE: setPoints can be within or outside of the actual 
%                         range of x, but must be monotonically increasing. Intended
%                         use is min(setPoints) = min(x), and max(setPoints) = max(x).
%                         If x contains values < min(setPoints) or > max(setPoints),
%                         those values will be colored corresponding to the low or high
%                         colors, respectively. 
%         'resolution'    number of points on the colormap returned. RGB values
%                         are not affected by this parameter. Equivalent to 
%                         e.g. the '100' in jet(100). Default: 100
%            
% Outputs: 
%
%     cinfo      structure containing:
%
%       rgb         numel(x) x 3 matrix of interpolated RGB colors (one for 
%                   each point in x).
%       cmap        resolution x 3 matrix of RGB values for custom colormap 
%                   (e.g. as returned by jet(100)). Each row defines the color 
%                   corresponding to the values in cmapVals (interpolated from
%                   x and setPoints, but not matching them). 
%       cmapVals    Values corresponding to the colors in cmap, interpolated from
%                   x and setPoints.
% 
% Usage:
%
% 
%     [cinfo] = customColormap(x,'resolution',100,'setPoints','colors',[-.2,median(x),0.05],colors = [0.568600000000000,0.0392000000000000,0.0392000000000000;1,1,1;0,0,0.568600000000000]);
%      rng = max([min(abs(vals)),max(abs(vals))]); [cinfo] = customColormap(x,'setPoints',[-rng,0,rng]);
% 
% MTJ 2020

    %% Default params
    
        res = 100;
        colors = [0,0,0.5686;...
                  1,1,1;...
                  0.5686,0.0392,0.0392]; % Blue (low vals), White (middle vals), Red (high vals)
        balance = 0;
        
    %% Parse params
    
        if ~isempty(varargin)   
            ind = find(strcmp('resolution',varargin),1);
            if ~isempty(ind)
                res = varargin{ind+1}; % reset resolution to passed val
            end
            ind = find(strcmp('colors',varargin),1);
            if ~isempty(ind)
                colors = varargin{ind+1}; % reset colors to passed val
            end
            ind = find(strcmp('setPoints',varargin),1);
            if ~isempty(ind)
                setPts = varargin{ind+1}; % reset setPoints to passed val
            end
            ind = find(strcmp('balance',varargin),1);
            if ~isempty(ind)
                balance = true;
            end
        end
    
    %% Interpolate colors
    
        % Interpolate evenly spaced map from n colors   
            
%             c1 = uisetcolor;
%             c2 = uisetcolor;
%             c3 = uisetcolor;           
        
        x = x(:);   % linearize
        
            % Interpolate
                numBaseColors = size(colors,1);
                
                mn = min(x);
                mx = max(x);
                
                % If setpoints are not provided, evenly space setpoints  
                % throughout range (min and max inclusive)
                    
                    if ~exist('setPts','var')
                        setPts = (mn: ( mx - mn ) / (numBaseColors-1) : mx)';                    
                    end

                % Now that we have setpoints, do we need them balanced?
                % Balance = endpoints expanded to be equidistant from median
                % If default setpoints are used, nothing should happen.
                
                    if balance
                        setPts(1) = median(setPts)-max( abs(setPts-median(setPts) ));
                        setPts(end) = median(setPts)+max( abs(setPts-median(setPts) ));
                    end
                    
                %% Interpolate RGB colormap                  
                
                    % Handle cases where scale is saturated on either or both ends:
                    
                        lSat = x < min(setPts);
                        uSat = x > max(setPts);
                        interpFor = ~or(lSat,uSat); % if there's no saturation, this has no effect
                    
                    % Generate the rgb vals for each point 
                    
                        rgb = zeros(numel(x),3);
                        rgb(interpFor,1) = interp1(setPts,colors(:,1),x(interpFor),'linear');
                        rgb(interpFor,2) = interp1(setPts,colors(:,2),x(interpFor),'linear');
                        rgb(interpFor,3) = interp1(setPts,colors(:,3),x(interpFor),'linear'); 
                    
                    % Generate the cmap and corresponding vals (in case it's needed)
                    
                        cmap = zeros(res,3);
                        cmapVals = (min(setPts): ( max(setPts) - min(setPts) ) / (res-1) : max(setPts))';
                        cmap(:,1) = interp1(setPts,colors(:,1),cmapVals,'linear');
                        cmap(:,2) = interp1(setPts,colors(:,2),cmapVals,'linear');
                        cmap(:,3) = interp1(setPts,colors(:,3),cmapVals,'linear'); 

                    % Add in saturations (if needed)
                    
                        if any(lSat)
                            rgb(lSat,:) = repmat(colors(1,:),sum(lSat),1);
                        end
                        if any(uSat)
                            rgb(uSat,:) = repmat(colors(end,:),sum(uSat),1);
                        end
                        
                    % Useful for checking our work:    
%                         [~,inds] = sort(x);
%                         rgb(inds,:);
                    
            cinfo.rgb = rgb;
            cinfo.cmap = cmap;
            cinfo.cmapVals = cmapVals;
            cinfo.setPts = setPts;

end

