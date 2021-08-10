function [newColors] = subsetColorStruct(colors,samplesubset)
%% Apply a spectrum subsetting operation to colors struct

%% 
    allinds = find(any(colors.inds_cat==samplesubset,2));
%     allinds = allinds(indscat);
    newColors.categories = colors.categories;
    newColors.inds_cat = colors.inds_cat(allinds);
    newColors.colorList = colors.colorList;
    newColors.rgb = colors.rgb(allinds,:);

end