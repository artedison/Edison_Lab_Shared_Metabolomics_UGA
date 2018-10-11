function [linInds,contours,windowInds,ppmtester] = getAdjacentTrajectories_simpler(matrix,currentppm,ridgeInds_row,ridgeInds_col,halfWindowWidth,viewWidth)
%% getAdjacentTrajectories_simpler

% Author: Michael T. Judge
% Version: 0.1
% Date: 04/12/2018

% Description:
%       Given peak position across time, this function builds a matrix of
%       indices that surround the peak. In other words, given a wandering
%       peak, take a wandering window around that peak. Accessory function
%       to ridgeCorrection.m
%
% Input: 
%       matrix: spectroscopic data, where spectra = rows, ppms = columns
%       currentppm: ppm vector for matrix
%       ridgeInds_row: row of matrix for each ridge point
%       ridgeInds_col: currentppm/matrix indices providing peak position in 
%           each row
%       halfWindowWidth: number of data points to include to left and right
%           of the ridge indicating peak position at different rows
%       viewWidth: spectral width considered for plotting (in ppms)
%
% Output: 
%       linInds: linear indices defining the peak position across spectra,
%           such that matrix(linInds) returns the ridge trajectory
%       contours: the result of matrix(linInds), giving matrix intensity at
%               each point in the window
%       windowInds: linear indices defining the window across spectra,
%           such that matrix(windowInds) returns the entire window
%       ppmtester: useful for extracting the ppm positions from windowInds,
%           accomplished by ppmtester(windowInds) 
%
% Log:
%
% Example run: 
%       See ridgeCorrection.m for context:
%       [linInds,contours,windowInds,ppmtester] = getAdjacentTrajectories_simpler(matrix,currentppm,ridgeInds_row,ridgeInds_col,halfWindowWidth,viewWidth)


%%
% Autoset the ROI:
    ROI = [min(currentppm(ridgeInds_col))-viewWidth,max(currentppm(ridgeInds_col)+viewWidth)];
    %ROI = [6,7];
    ROIinds = matchPPMs(ROI(1),currentppm):matchPPMs(ROI(2),currentppm);
            % Get the linear indices of the 
                linInds = sub2ind(size(matrix),ridgeInds_row',ridgeInds_col);
%             % Check to see that these are right
%                 figure, hold on
%                     surf(currentppm(ROIinds),1:size(matrix,1),matrix(:,ROIinds),'FaceColor','Interp');%,hold on
%                     shading interp
%                     xlabel('trajectory')
%                     zlabel('Scaled Intensity')
%                     ylabel('Time (h)')      
%                     set(gca,'xdir','reverse')
%                     scatter3(currentppm(ridgeInds_col), ridgeInds_row, matrix(linInds))
%                 hold off 
            % Get the linear indices of the lines (adjacent trajectories)
                % Issue: ridges are
                %           - discontinuous
                %           - different lengths
%                 for i = 1:length(min(ridgeInds_row):max(ridgeInds_row))
%                     
%                 end;.;
                
                ladjusts = fliplr(1:halfWindowWidth).*size(matrix,1); 
                radjusts = (1:halfWindowWidth).*size(matrix,1);
                
%                 lcontours = repmat(ridgeInds_col',1,halfWindowWidth)-repmat(ladjusts,size(ridgeInds_col',1),1);
%                 rcontours = repmat(ridgeInds_col',1,halfWindowWidth)+repmat(radjusts,size(ridgeInds_col',1),1);
                lcontours = repmat(linInds',1,halfWindowWidth)-repmat(ladjusts,length(linInds'),1);
                rcontours = repmat(linInds',1,halfWindowWidth)+repmat(radjusts,length(linInds'),1);
                
                windowInds = [lcontours,  linInds', rcontours];    % these are linear indices
                            % See if this is what we want 
                                 ppmtester = repmat(currentppm,size(matrix,1),1);
%                                 figure,hold on,plot(ppmtester(windowInds),1:length(fridge),'k','linewidth',0.5)
%                                 scatter(ppmtester(ridgeInds),1:length(fridge),'r')
            % Use the linear index for the entire matrix to pull out the
            % intensities along the trajectory and adjacent trajectories:
                contours = matrix(windowInds);
            
            % Visualize these trajectories as a straightened (mapped) peak and 
%             %   including the surrounding regions:                    
%                 figure, hold on
%                         surf(1:size(contours,2),1:size(contours,1),contours,'FaceColor','Interp');%,hold on
%                         %shading interp
%                         xlabel('trajectory')
%                         zlabel('Scaled Intensity')
%                         ylabel('Time (h)')
%                         title(['Mapped trajectory window for Trajectory ',num2str(i)])
%                         %figuresaver('fig')           
%                         %plot(ppm(windowInds),'k','linewidth',1) 
%                         plot3(repmat(halfWindowWidth+1,size(contours,1),1),timepoints,rvals,'r','linewidth',3)
%                         %scatter3(1:length(fridge),ppmtester(ridgeInds),countours,'r')
%                         hold off

%%  
%             % Check to see that these bounds are right
%                 figure, hold on
%                     surf(currentppm(ROIinds),1:size(matrix,1),matrix(:,ROIinds),'FaceColor','Interp');%,hold on
%                     shading interp
%                     xlabel('trajectory')
%                     zlabel('Scaled Intensity')
%                     ylabel('Time (h)')
%                     title('Uncorrected Ridge on Unsmoothed Data')
%                     set(gca,'xdir','reverse')
%                     % Plot the ridge in blue:
%                       %scatter3(currentppm(ridgeInds_col), ridgeInds_row, matrix(linInds))
%                       plot3(currentppm(ridgeInds_col), ridgeInds_row, matrix(linInds),'b','linewidth',2)
%                     % Plot the bounds in red:
%                       [endRows,~] = ind2sub(size(matrix),[windowInds(:,1),windowInds(:,end)]);
%                       plot3(ppmtester(windowInds(:,1)),endRows(:,1),matrix(windowInds(:,1)),'r','linewidth',3)
%                       plot3(ppmtester(windowInds(:,end)),endRows(:,2),matrix(windowInds(:,end)),'r','linewidth',3)                    
%                 hold off 
                        
end