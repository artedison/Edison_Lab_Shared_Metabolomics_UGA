function [integral] =integrate_AJB(X_in,Y_in,ppm_in,cursor_info,mut_ID_num)
%
% Integrate the area under the NMR curve between designated right and left ppm boundaries for all X-value arrays in X_input, subtracting the baseline area (line draw between two boundries) from total area. 
%
% Arguments:
% 
% X_in            N x P matrix of spectral feature intensities (P) for each sample (N)
%
% Y_in            Y_vector corresponding to Sample Identities
%
% ppm_in         They array containing all ppm input values
% 
% cursor_info    The cursor info table generated after exporting selected
%                points between which you wish to integrate
%                 
% mut_ID_num     The ID number for the data set designated as mutant 
%
% Output:
%
% [Integral]
%      > A structure Field called summary with an array listing:
%           -Peak identities
%           -WT_mean
%           -WT_Stdev
%           -Mut_mean
%           -Mut_Stdev
%           -Fold-change 
%           -p-value
%           -correct p-value (q-value) as obtained from the Benjamini-Hochberg
%               method (1995)

%      > Fields using the highest peak ppm as identifiers with the following: 
%                   - WT_area: Mean integration area value for the wild-type
%                             (reference) dataset
%                   - WT_StDev: Standard deviation value for the wild-type
%                             (reference) dataset
%                   - Mutant_area: Mean integration area value for the
%                                 mutant(test) dataset
%                   - Mutant_StDev: Standard deviation value for the mutant
%                                  (test) dataset
%                   - Fold_Change: (Mutant_area/WT_area) or (test/control)--
%                                 this is a 1-centered fold change value
%                   - P_Value: Uses Student's two-sample t-Test (ttest2) to
%                             determine a p-value for whether WT_area and 
%                             Mutant_area significantly differ from eachother
%                   - Data_Table: Place the values determined above into an
%                                   array for easy downstream applications in the order
%                                   (WT_area; WT_StDev; Mutant_area; Mutant_StDev;
%                                    Fold_Change; P_Value
%                   - Raw_Info: A structure array containing the data used
%                               to determine the values above.
%                              Structure Fields:
%                                  * Interval: The initial ppm values selected by the 
%                                               user to define peak bounds used in the program
%                                  * WT_replicates: Peak area determined for each wild-type 
%                                                   (control) replicate
%                                  * Mutant_replicates: Peak area determined for each mutant 
%                                                   (test) replicate
%                                  * MT_ppm_Bounds: Specific peak bounds (ppm) determined by the
%                                                   program for each mutant (test) replicate
%                                  * WT_ppm_Bounds: Specific peak bounds (ppm) determined by the
%                                                   program for each wild-type (control) replicate
%                                  * WT_Baseline_Slope: Slope of the baseline determined between 
%                                                       the wild-type(control) bounds for each replicate
%                                  * MT_Baseline_Slope: Slope of the baseline determined between 
%                                                       the mutant (test)bounds for each replicate
%      
%            Also calls "Plot_Integral_AJB" which creates an "Integrals" Folder in current directory that
%            contains figures highlighting the peak integrated and plotting the baseline used for peak area 
%            determination for every wild-type (control) and mutant (test)replicate. Total area and baseline 
%            area are displayed on each figure 
%
%
% Make arrays from the "cursor_info" data table for the Left and Right Bounds
% of the peak to be integrated
%
[ppmLeft,ppmRight] = PeakBounds(cursor_info);
%
% Find out how many samples are included in the X_input matrix
[rowX,colX]=size(X_in);
%
% Find unique values and their indices for the Y-vector (which lists Sample Identities)
%
%       This finds the unique values within 'Y_in' and the index vectors 'ia'
%       and 'ic' such that C = Y_in(ia) and Y_in(ic) = C
%       Example: 
%          A = [9 2 9 5];
%            [C, ia, ic] = unique(A)
%          C = 1×3    [2,5,9]
%          ia = 3×1   [2;4;1]
%          ic = 4×1   [3;1;3;2]
[C,ia,ic]=unique(Y_in);
%
% How many unique sample identities exist? (should only be 2 for wild-type
%                                           vs. mutant)
[uniqueY,~]=size(C);


%   Go through the below process for every peak selected
for k=1:length(ppmLeft)
    
%   Must flip ppmLeft and ppmRight if ppmLeft is bigger than ppmRight
    if ppmLeft(k)>ppmRight(k)
        ppmLflip = ppmRight(k);
        ppmRflip = ppmLeft(k);
        ppmLeft(k) = ppmLflip;
        ppmRight(k) = ppmRflip;
        clear ppmLflip;
        clear ppmRflip;
    end
    xvect=[];
    
%   Find the column numbers corresponding to the left and right bounds
    [~,idxL]=min(abs(ppmLeft(k)-ppm_in));
    [~,idxR]=min(abs(ppmRight(k)-ppm_in));
   
%   Pull out only the intensity values for the ppms between left and right
%       bounds
        Xseg=X_in(:,idxL:idxR);
        
%   Find the highest peak within the middle of the range to act as an identifier
        for z = idxL:idxR;
            xvect=[xvect, mean(X_in(1:rowX,z))];
        end
        [~,sz_xvect]=size(xvect);
        x_vect_fth = round(sz_xvect/5);
        peakVal = max(xvect(x_vect_fth:x_vect_fth*4));
        peakppm(k)=ppm_in(find(xvect==peakVal)+idxL);
        [~,idxP]=min(abs(peakVal-Xseg(k,:)));
        
% Go through the below process for every sample
    for pl = 1:rowX
       
%   Find the highest peak in the middle of the range for each specific sample
        sam_P = max(Xseg(pl,x_vect_fth:x_vect_fth*4));
        [~,idxSam_P]=min(abs(sam_P-Xseg(pl,:)));
        
%   Find the two valleys that bound the peak of interest for each individual sample 
%       (look in from each boundry) to account for variation between samples
        Left_side = Xseg(pl,1:idxSam_P);
        Lmin = min(Left_side);
        Right_side = Xseg(pl,idxSam_P:end);
        Rmin = min(Right_side);
        new_idxL=find(Xseg(pl,:)==Lmin);
        new_idxR=find(Xseg(pl,:)==Rmin);
        
%   If Multiple Points are LB then trim down to the first point
                    if size(new_idxL,2)>1
                        new_idxL = new_idxL(1);
                    end
        
%   If Multiple Points are RB then trim this down to the last point
                    if size(new_idxR,2)>1
                        new_idxR = new_idxR(end);
                    end
                    
%   Using the peak bounds specific to this particular sample, create new 
%       peak intesity and ppm arrays                    
        new_Xseg = Xseg(pl,new_idxL:new_idxR);
        LB = idxL+(new_idxL-1);
        RB = idxR-(size(Xseg,2)-new_idxR);
        new_ppm = ppm_in(LB:1:RB);
    
%   Create a peak intesity array for the baseline that connects the two peak bounds
%       by a straight line
    Y1 = new_Xseg(1);
    Y2 = new_Xseg(end);
    Xdiff = size(new_Xseg,2);
    Ydiff = Y2 - Y1;
    Slope_Baseline(pl,k) = (Ydiff)/(Xdiff);
    clear BaseLine_Xseg;
    for lng = 1:size(new_Xseg,2)
        BaseLine_Xseg(lng) = new_Xseg(1)+((lng-1)*Slope_Baseline(pl,k));   
    end
    
%   If the slope at the ends of the baseline increase/decrease faster/slower than 
%       the slope for the actual peak data, trim the peak boundries used by 1 
        if BaseLine_Xseg(2)>new_Xseg(2) || BaseLine_Xseg(end-1)>new_Xseg(end-1) && size(new_Xseg,2)>3
            
%           Trim left bound if baseline is higher than peak line
                if BaseLine_Xseg(2)>new_Xseg(2)
                    new_Xseg = new_Xseg(2:end);
                    new_ppm = new_ppm(2:end);
                end
                
%           Trim right bound if baseline is higher than peak line
                if BaseLine_Xseg(end-1)>new_Xseg(end-1)
                    new_Xseg = new_Xseg(1:end-1);
                    new_ppm = new_ppm(1:end-1);
                end
                
%           Obtain a new baseline array using the trimmed bounds
                Y1 = new_Xseg(1);
                Y2 = new_Xseg(end);
                Xdiff = size(new_Xseg,2);
                Ydiff = Y2 - Y1;
                Slope_Baseline(pl,k) = (Ydiff)/(Xdiff);
                clear BaseLine_Xseg;
                    for lng = 1:size(new_Xseg,2)
                        BaseLine_Xseg(lng) = new_Xseg(1)+((lng-1)*Slope_Baseline(pl,k));   
                    end
        end
        
%   After setting peak bounds, If X_seg isn't over baseline for more than 4 pts, then the 
%       peak isn't detectable above baseline print error message
        if size(new_Xseg,2)<=4 
            Left_Bounds(pl,k) = 0;
            Right_Bounds(pl,k) = 0;
            Slope_Baseline(pl,k) = 0;
            baseline_area(pl,k) = 0;
            total_area(pl,k) = 0;
            integralX(pl,k) = 0;
            SMPL = num2str(pl);
            Peak1 = num2str(ppmLeft(k));
            Peak2 = num2str(ppmRight(k));
            ermsg = strcat('The peak between -', Peak1,'- and -', Peak2, '- for sample -',SMPL,'- cannot be seen above baseline \n \n');
            fprintf(ermsg);
            continue
        end
        
%   If baseline value is ever greater than data value after assigning the peak bounds, baseline value = data value
        for BL_Check = 1:size(new_Xseg,2)
                if BaseLine_Xseg(BL_Check) > new_Xseg(BL_Check)
                        BaseLine_Xseg(BL_Check) = new_Xseg(BL_Check);
                end
        end
        
%   Catalog the left and right boundries used for every sample after trimming
        Left_Bounds(pl,k)= new_ppm(1);
        Right_Bounds(pl,k)= new_ppm(end);
        
%   Integrate baseline area(trapz method) between two valley bounds for the current sample
    baseline_area(pl,k)=trapz(BaseLine_Xseg); 
    
%   Integrate total area(trapz method) between two valley bounds for the current sample
    total_area(pl,k)=trapz(new_Xseg);
    
%   Find the peak area between two valley points by (Total area - Baseline Area) 
    integralX(pl,k)=total_area(pl,k)-baseline_area(pl,k);
    end 
end

%   Make a structure table containing 2 fields (wild-type and mutant) and populate 
%       these fields using a matrix with rows representing all the samples for that field
%       and columns representing the different peaks for which integrals were calculated  
        for l=1:uniqueY
            ident = num2str(C(l));
            S = 'Sample';
            I = 'integral';
            Sident = strcat(I,'_',S,'_',ident);
            Integral_Table.(Sident)= integralX(Y_in == C(l),:);
        end
        
        Integral_Table_cell = struct2cell(Integral_Table);
        
%   Figure out if mutants IDs are provided first or second and set up following scripts accordingly  
        if mut_ID_num == C(1)
            integral_MT = Integral_Table_cell{1};
            integral_WT = Integral_Table_cell{2};
        else
            integral_MT = Integral_Table_cell{2};
            integral_WT = Integral_Table_cell{1};
        end
        
%   Start building the integral table
        [row,col]=size(integralX);
        ave_integral_MT=zeros(1,col);
        std_MT=zeros(1,col);
        ave_integral_WT=zeros(1,col);
        std_WT=zeros(1,col);
        FC = zeros(1,col);
        pval = zeros(1,col);
        for l = 1:col
            ave_integral_MT(l) = mean(integral_MT(:,l));
            std_MT(l) = std(integral_MT(:,l));
            ave_integral_WT(l) = mean(integral_WT(:,l));
            std_WT(l) = std(integral_WT(:,l));
            FC(l) = ave_integral_MT(l)/ave_integral_WT(l);
        end
        for p = 1:col   
            [H,P] = ttest2(integral_MT(:,p),integral_WT(:,p));
            pval(p) = P;
            
        end
%
%   Find the corresponding Benjamini-Hochberg (1995) corrected p-values (q-values) for
%   all the samples (this calls a function from the bioinformatics toolbox)
       
            qval = mafdr(pval, 'BHFDR',true);
            Tab_labels = ["Peak_ppm", "WT_mean", "WT_Stdev", "Mut_mean", "Mut_Stdev","Fold_Change", "p_value", "FDR_value"];
            stats_table = table([peakppm;ave_integral_WT;std_WT;ave_integral_MT;std_MT;FC;pval;qval],'RowNames',Tab_labels);

            
       


%   Finish by finding out how many different peaks were integrated and then fill in the
%   outputs and generate the figures peak by peak
      integral.Summary_Table=stats_table;
        [ppmX, ppmY]=size(peakppm);
        for l =1:ppmY
            Raw_Info.Interval = [ppmLeft(l),ppmRight(l)];
            Raw_Info.WT_replicates = integral_WT(:,l);
            Raw_Info.Mutant_replicates = integral_MT(:,l);
            Raw_Info.MT_ppm_Bounds = [Left_Bounds(Y_in == mut_ID_num ,l),Right_Bounds(Y_in == mut_ID_num,l)];
            Raw_Info.WT_ppm_Bounds = [Left_Bounds(Y_in ~= mut_ID_num ,l),Right_Bounds(Y_in ~= mut_ID_num,l)];
            Raw_Info.WT_Baseline_Slope = Slope_Baseline(Y_in ~= mut_ID_num,l);
            Raw_Info.MT_Baseline_Slope = Slope_Baseline(Y_in == mut_ID_num,l);
            Info.WT_area = ave_integral_WT(l);
            Info.WT_StDev = std_WT(l);
            Info.Mutant_area = ave_integral_MT(l);
            Info.Mutant_StDev = std_MT(l);
            Info.Fold_Change = FC(l);
            Info.P_Value = pval(l); 
            Info.Data_Table = [ave_integral_WT(l);std_WT(l);ave_integral_MT(l);std_MT(l);FC(l);pval(l)];
            Info.Raw_Info = Raw_Info;
            field = num2str(peakppm(l));
            field = strrep(field,'.','_');
            P = 'Peak';
            Sfield = strcat(P,'_',field);
            integral.(Sfield) = Info;
            

%   Call the "Plot_Integral_AJB" Program to generate figures zooming in on
%       the peak used for integration and plotting the baseline for every sample
%       and every peak integrated
            Plot_Integral_AJB(integral,Sfield,X_in,ppm_in,Y_in,mut_ID_num);
        end


