function [itpeaks] = Peakpick1D_per_spectra(matrix,ppm,peakthresh,mode, representative_spectra)
    
%% Per sample peak picking 
% uses PeakPick1D_noFigures to supress figure output

% Nofig function call
% [peakmatrix,shifts]=Peakpick1D_noFigures(X,ppm,represent,peakthresh,mode)

% Original function
%[peakmatrix,shifts]=Peakpick1D(X,ppm,represent,peakthresh,mode)

% This functions creates a loop based on the dataset to carryout peak picking per spectra
% Does not need a representative spectra therefore the parameter "represent" is set to 1
%
% Input: 
%           matrix           -  where rows are samples and columns are ppms
%           ppm              -  ppm vector with the same number of columns as matrix and
%                                    one row
%           peakthresh  - This is data dependent - should be low enough to
%                                    peak as many peaks as possible - use
%                                    Peakpick1D to optimize threshold
%                                    beforehand
%           mode            - see Peakpick1D
%           representative_spectra          - Spectra user wants to see plotted (for diplaying only)
%                                                               input is a number (double)
%                                                               indicating the index of the spectra
%                                                               of choice - if none selected then
%                                                               the middle spectra of the data set
%                                                               is displayed

% Output:            
%           itpeaks          a structure containing as fields:
%                                  - ints - intensity for each spectra
%                                  - shifts - chemichal shift (ppm) vector
%                                  for each spectra
%
% Note: Structures names can be user defined, however substructures and
% fields need to remain as defined by each function 

% example parameters setup:
%              matrix = wrk_data.XRBA;
%         ppm = wrk_data.ppmR;
%         peakthresh = 0.2
%         mode = 'Complex'


             itpeaks = struct();
        
             for i=1:size(matrix)
                 
             [itpeaks.iteration(i).ints, itpeaks.iteration(i).shifts]= Peakpick1D_noFigures( matrix(i,:),ppm,1,peakthresh,mode); %interger 1 for a single spectra
      
             end

             
%% displaying a representative spectra 

            if nargin < 5 || isempty(representative_spectra)
                
                  rsp = round(size(matrix,1)/2); %random spectra
                  
            else
                     rsp = representative_spectra;
             end
             
               figure, 
             hold
             plotr ( ppm, matrix ( rsp ,:) , 'k')
             plotr( itpeaks.iteration(rsp).shifts,...
                        itpeaks.iteration(rsp).ints, 'ob') 
             
             title([ 'Displaying spectra n: ', string(rsp)])
         
             
             
             