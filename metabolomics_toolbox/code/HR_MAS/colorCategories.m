function [colors] = colorCategories(colorBy,varargin)

% 'cmap'      'jet', 'lines', etc
% 'bgColors'  {'w','k'} or [1,1,1;0,0,0]
    

    %% Parse inputs
    
        % background color(s) - colors to avoid (passthrough to
        % distinguishable_colors()
    % defaults
        bgColors = [1 1 1; 0 0 0];
        
    if ~isempty(varargin)
        ind = strcmp('cmap',varargin);
        if any(ind)
            cmap = varargin{find(ind,1)+1};
        end
        
        ind = strcmp('bgColors',varargin);
        if any(ind)
            bgColors = varargin{find(ind,1)+1};
        end
        
    end
        
    %% Interpolate colors

    
        % If x is categorical (can include multiple columns)
                            
            colors = struct(); % We'll keep the output in a struct object

            % If we want to put the categories on a continuous scale:
        
                % Not sure what to do here yet
                
            % If we want distinguishable colors:
            
                % If colorBy is not a table, cast to table
            
                if ~istable(colorBy)
                    colorBy = table(colorBy);
                end
            
            [colors.categories,~, colors.inds_cat] = unique(colorBy,'rows','stable');
            if exist('cmap','var')
                eval(sprintf('colors.colorList = %s(height(colors.categories))',cmap));
            else
                colors.colorList = distinguishable_colors(height(colors.categories),bgColors);
            end
            colors.rgb = colors.colorList(colors.inds_cat,:);
            
end

