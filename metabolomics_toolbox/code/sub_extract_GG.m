function [new_X, new_metadata] = sub_extract_GG( X, metadata, metadata_field, metadata_element, remove);

  % Author: Goncalo Gouveia
    % Version: 0.1
    % Tested on Matlab Version R2019a
    % Date: 8Jun2020
    % 
    %% Description:
    %       Given that the matrix X (rows are samples and columns are
    %       variables) is row matched to the metadata file (every row in
    %       the metadata corresponds to a row in the matrix X),
    %       sub_extract_GG creates a new matrix and corresponding metadata
    %       matrix based on the metadata element specified
    %
    % 
    % 
    %% Input:
    %       X        : 
    %       Numerical matrix where rows are samples and columns
    %       are variables
    %       Also works with structures if following the same order (each
    %       row is a sample )
    %
    %%       metadata     :  
    %      table where rows are samples (matched to X) and
    %       columns are metadata variables
    %       Can be a single vector too, but MUST BE A STR OR CHR
    %
    %%       metadata_field   :  
    %             the column field of the metadata selected
    %             to subsample the data from.
    %             Can be both a vector or a column in a matrix or table 
    %             i.e. metadata.column_to_select_name
    %
    %
    %%       metadata_element   :  
    %       this needs to be a str or character.  The
    %       function will match partial matches. 
    %       example: if the user only
    %       wants samples of the fruit type "apple" in the metadata_field containing
    %       [ 'apple','banana','orange', 'pineapple' etc] if the metadata_element 'apple'
    %       is inputed, it matches to all rows containing 'apple' and 'pineapple'
    %
    %%       remove   :  
    %       If set to "remove" then the rows matching the metadata_element
    %       will be removed 
    %       If set to 'keep' the rows
    %       matching the metadata_element will be extracted to a new
    %       variable
    %       
    %
    %% Info:
    %       Does not work if the metadata table field is a number/double
    %       If there are elements present that are similar it will also
    %       select that, example: search for "bana" and it will find
    %       "banana" but also "abanana" "bananas" "banana1" "banana100" etc
    %       
    %       when calling the function suppress the output by using  ->   "  ;  "
    %       at the end of the function call
    %       i.e. [wrk_data.Xnew, Tnew] = sub_extract_GG ( wrk_data.X, Tmeta, Tmeta.fraction_bio, 'b_', 'keep') ;
    
    
metadata_field = string( metadata_field);
    
if ismember( string('remove'), remove) == 1
    
new_X = X( ~contains(metadata_field,metadata_element),:);

new_metadata = metadata( ~contains(metadata_field,metadata_element),:);


echo off;
warning ([ 'removing rows containing ->  ', metadata_element])

else
    
 new_X = X( contains(metadata_field,metadata_element),:);

new_metadata = metadata( contains(metadata_field,metadata_element),:);   

echo off;
metadata_element = char(metadata_element);

warning ([ 'extracting rows containing ->  ', metadata_element])

end


end




