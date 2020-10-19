function [updatedStruct] = removeRows(dataStruct,expectedLength,rowsToRemove)
%% removeRows

    % Author: MTJ
    % Version: 0.1
    % Tested on Matlab Version R2020a
    % Date: JUN2020
    %
    % Description:
    %       
    %       Removes rows (e.g. outlier samples) in structure elements with 
    %       the same number of rows as the X matrix (following a 
    %       structure data storage scheme). Keeps metadata and spectra in 
    %       concordance instead of keeping track of indices. Searches
    %       contents of dataStruct for matrices or tables with number of
    %       rows = expectedLength, then removes rows following indices 
    %       provided in rowsToRemove. For example, it accomplishes the
    %       following in a single command:
    %           dataStruct.X(rowsToRemove,:) = [];
    %           dataStruct.XR(rowsToRemove,:) = [];
    %           dataStruct.XLSN(rowsToRemove,:) = [];
    %           dataStruct.metadata(rowsToRemove,:) = [];
    %       The function ignores fields whose contents do not have
    %       expectedLength rows. For example: 
    %           If size(dataStruct.ppm) is [1,30000], and expectedLength =
    %           42, then the function will simply pass through dataStruct.ppm
    %           without modification. 
    %
    % Input:
    %
    %       dataStruct:     Structure containing spectra and related data, e.g. fields:
    %                           X:          Spectral matrix            e.g. size(X)  = [42,30000]                    
    %                           XR:         X, regions removed         e.g. size(XR) = [42,25000]    
    %                           XRSN:       X, scaled and normalized   e.g. size(XR) = [42,25000]    
    %                           metadata:   table of metadata, rows correspond with rows
    %                                       of X                       e.g. height(metadata) = 42                                    
    %                           ppm:        ppm vector for X (ignored) e.g. size(ppm)   = [1,30000]    
    %                           
    %                           ppmR:       ppm vector for XR (ignored) e.g. size(ppmR) = [1,25000]
    %
    %       expectedLength: number of rows in the variable to be indexed.
    %                       Typically size(X,1) or height(metadata)
    %
    %       rowsToRemove:   indices of rows to delete. Can be logical
    %                       vector or array of integers. Typically these
    %                       are indices of outliers. 
    %
    % Output:
    %       
    %       updatedStruct:  returned version of dataStruct with rows
    %                       removed for appropriate fields. Also contains a
    %                       new substruct field:
    %               removeRows.expectedLength:  passthrough, record of param
    %               removeRows.rowsToRemove:    passthrough, record of param
    %                   NOTE: if removeRows was applied multiple times, the
    %                   indices will be concatenated to preserve the
    %                   previous record
    %
    % Log:
    %
    %
    %
    % Example run:
    %
    %       [projectpr1d_OR] = removeRows(projectpr1d,height(projectpr1d.metadata),projectpr1d.outlierInds);
    %   
 
%% Find the fields with rows 

     updatedStruct = dataStruct;
%     expectedLength = height(projectpr1d.metadata);
%     rowsToRemove = logical(projectpr1d.metadata.outliers);
%     [updatedStruct] = removeRows(updatedStruct,expectedLength,rowsToRemove)

    if islogical(rowsToRemove)
        rowsToRemove = find(rowsToRemove);
    end
       fields = fieldnames(updatedStruct);  
       
    % Don't include previous inds in this (in case of 62 rows)
        if isfield(updatedStruct,'removedRows')
            fields(strcmp(fields,'removedRows')) = []; 
        end
        
    % Remove the rows in the new struct if their variable are the right
    % length (e.g. if they're actually indexed by rowsToRemove)
       for i = 1:length(fields)
           if size(updatedStruct.(fields{i}),1) == expectedLength
                updatedStruct.(fields{i})(rowsToRemove,:) = [];
           else
                fprintf(['\n\n\tremove_rows: the field ',fields{i},' was not modified.\n'])
           end
       end
       
    % Store the params we used in the returned struct
    if isfield(updatedStruct,'removeRows')
        updatedStruct.removeRows.rowsToRemove = [updatedStruct.removedRows;rowsToRemove(:)];
        updatedStruct.removeRows.expectedLength = expectedLength;
    else 
        updatedStruct.removeRows.rowsToRemove = rowsToRemove(:);
        updatedStruct.removeRows.expectedLength = expectedLength;
    end
    
end