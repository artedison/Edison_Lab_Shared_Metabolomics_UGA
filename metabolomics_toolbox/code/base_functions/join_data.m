function [newX,newPpm] = join_data(X,ppm,varargin)

% Author: Rahil Taujale
% Version: 0.2
% Date: 06/26/2017

% Description:
%       Simulates a Gaussian curve for each data point so that each data 
%       point (each column) ends up being 50 columns that have a Gaussian peak
%       Height of the peak corresponds to the intensity of the data point.
%       Intensity is scaled using the max intensitites of the matrices.
%       Adds this dataset to X along with a corresponding ppm axis.
%
% Input: 
%       X       : stack of 1D spectra
%       ppm     : chemical shift vector
%       varargin: multiple data matrices that you want to add to X in form
%       of a Gaussian peak. Number of rows in each of these matrices must
%       match the no. of rows in X.
%
% Output: 
%       A new matrix with added data points.
%       A new ppm axis where simulated ppm values are added for every new 
%       data point added. 
%
% Log:
%       Ver 0.2 - changed peak widht spacing from 50 to 10 for simulated
%       Gaussian peaks.
%
% Example run: [JointX,JointPpm]=join_data_new(XALN,ppmR,new_peaks);

    x=-9:1:10;
    y=gaussmf(x,[3,0]);
    newX=X;
    [A,B]=size(X);
    newPpm=ppm;
    for n=1:length(varargin)
        data_sc=scale_cross_platform(X,varargin{n});
        [samples,peaks]=size(varargin{n});
        if (A~=samples)
            error('Number of samples(rows) do not match across all matrices.')
        end
        for j = 1:peaks
            for i = 1:samples
                Add_peak(i,:,j)= data_sc(i,j)*y;
            end
            newX = [newX Add_peak(:,1:end,j)];
            ppmgap = (newPpm(1,end)-newPpm(1,end-1));
            [a,b,c] = size(Add_peak);
            Add_ppm = zeros(1,b);
            for k= 1:b
                Add_ppm(1,k)= newPpm(1,end)+ppmgap*k;
            end       
            newPpm=[newPpm Add_ppm];
        end
    end
end