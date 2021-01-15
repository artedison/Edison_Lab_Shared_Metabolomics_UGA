function [clrGrps, colorCats] = colorData(colorBy)

    %% Parse inputs
        % *** need to pass customColormap vals through this function
        
        % *** need to parse out whether continuous or categorical
        %x = aucs_proms;
    
    %% Interpolate colors
    
        % Interpolate evenly spaced map from n colors   
            
%             c1 = uisetcolor;
%             c2 = uisetcolor;
%             c3 = uisetcolor;           
            % If x is continuous:
                            
          %      [rgb, cmap, cmapVals] = customColormap(colorBy,'colors',colors);

                % Can we just pass a native MATLAB colormap as colors above, or must do the following?
                
                    % % Fit to colormap
                    %     % Scale values
                    % 
                    %         xs = round((x- min(x(:)) )  /...
                    %             ( max(x(:) )-min( x(:) ) ),2);
                    %         xs = int8(xs*100+1);
                    % 
                    %     cmap = jet(101); 
                    % 
                    %     rgb = cmap(xs,:);
                    %     vals = [];
            
        % By row? 
        
            % Just normalize the rows first
        
        % If x is categorical
        
%             [clrGrps, colorCats] = unique(colorBy,'rows');
%             distinguishable_colors(clrGrps)
            [colors] = colorCategories(colorby);
            clrGrps = colors;
            colorCats = colors;
        
end

