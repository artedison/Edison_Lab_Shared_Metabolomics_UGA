function [colors] = colorCategories(colorBy)

    

    %% Parse inputs
    
        % background color(s) - colors to avoid (passthrough to
        % distinguishable_colors()
    
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
            
            [colors.categories,~, colors.inds_cat] = unique(colorBy,'rows');
            colors.colorList = distinguishable_colors(height(colors.categories),[1 1 1; 0 0 0]);
            colors.rgb = colors.colorList(colors.inds_cat,:);
            
end

