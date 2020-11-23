function [] = highlightROIs(ROIs,height,varargin)

%% highlightROIs
%{
    This function highlights the ppm regions provided in 'ROIs'.
        % Note: ROIs must be 2 x n 
    Hope this helps,
    MJ 2MAY2017
%}

edgeColor = 'none';
color = 'r';
transparency = 0.1;
%colorVect = {color};      % in the absence of being set explicitly, this is a passthrough

    % Read in any options:
        if ~isempty(varargin)
            varnames = varargin(1:2:end);
            varvals = varargin(2:2:end);

            varInd = strcmp(varnames,'color');
            if any(varInd)
                color = varvals{varInd};
                colorVect = repmat({color},size(ROIs,2),1);
            end
            
            varInd = strcmp(varnames,'transparency');
            if any(varInd)
                transparency = varvals{varInd};
            end
            
            varInd = strcmp(varnames,'edgeColor');

            if any(varInd)
                edgeColor = varvals{varInd};
            end                
           
            varInd = strcmp(varnames,'colorVect');
            
            if any(varInd)
                colorVect = varvals{varInd};
                colorVect = num2cell(colorVect,2);
            end                
            
            
            % StackSpectra params:
            varInd = contains(varnames,'horzshift');
            if any(varInd)
                stackParams.horzshift = varvals{varInd};
            end
            varInd = contains(varnames,'vertshift');
            if any(varInd)
                stackParams.vertshift = varvals{varInd};
            end    
            varInd = contains(varnames,'numberOfSpectra');
            if any(varInd)
                stackParams.numberOfSpectra = varvals{varInd};
            end    
            varInd = contains(varnames,'extension');
            if any(varInd)
                stackParams.numberOfSpectra = varvals{varInd};
            end    
        end
        
        if ~exist('colorVect','var')
            colorVect = repmat({color},size(ROIs,2),1); % in the absence of being set explicitly above, this is a passthrough
        end
        
    %% Make the highlights
        for i = 1:size(ROIs,2)
            %% Calculate the coordinates
                width = (ROIs(2,i)-ROIs(1,i));
                llx = ROIs(1,i); % 
                lux = llx;
                rlx = llx + width; % = rux
                rux = llx + width; % = rux;
                lly = 0;
                luy = height;
                rly = 0;
                ruy = height;
                xcoords = [llx rlx rux lux];
                ycoords = [lly rly ruy luy];
                
            %% If it's a stackSpectra plot, modify the coordinates:
                if exist('stackParams','var')
                    
                    % luy and ruy coords need to shift based on vertshift, numberOfSpectra
                        maxVertShift = stackParams.vertshift * stackParams.numberOfSpectra;
                        luy = 0 + stackParams.vertshift * stackParams.extension;
                        ruy = 0 + stackParams.vertshift * stackParams.extension;
                        lly = lly - maxVertShift - stackParams.vertshift * stackParams.extension;
                        rly = rly - maxVertShift - stackParams.vertshift * stackParams.extension;
                        ycoords = [lly rly ruy luy]; 
                        
                    % lux and rux coords need to shift based on horzshift, numberOfSpectra
                        maxHorzShift = stackParams.horzshift * stackParams.numberOfSpectra;
                        llx = llx + stackParams.horzshift * stackParams.extension;
                        rlx = rlx + stackParams.horzshift * stackParams.extension;
                        lux = lux - maxHorzShift;
                        rux = rux - maxHorzShift;
                        xcoords = [llx rlx rux lux]; 
                end

            %% Make the patch:
                p=patch(xcoords,ycoords,'r');

            %% Modify the patch in various ways:
                % Face Color
                    
                    set(p,'FaceColor',colorVect{i});

                % Edge Color
                
                    set(p,'EdgeColor',edgeColor);
                    
                % Transparency

                    set(p,'FaceAlpha',transparency);

        end
end