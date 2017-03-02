function [drivers] = interactions2drivers(interactions)

% 
drivers = unique( sort( reshape(interactions,size(interactions,1),[]) ) );
% matrix -> targetFile_secondary 

end