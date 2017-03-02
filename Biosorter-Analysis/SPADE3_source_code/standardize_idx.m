function idx_new = standardize_idx(idx)

idx_new = zeros(size(idx));
possible_idx = setdiff(unique(idx),0);
for i=1:length(possible_idx)
    idx_new(idx==possible_idx(i)) = i;
end


