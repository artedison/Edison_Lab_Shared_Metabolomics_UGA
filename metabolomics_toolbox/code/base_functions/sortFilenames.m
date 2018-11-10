function filenames_new=sortFilenames(filenames)
% Sort file names by natural ordering using numbers in the filename.
% Sort by leading digits, separated from text by an underscore
% For example filenames = {'1_a','11_a','2_a'}
% If encounters an error, return the input unchanged
% filenames: input file names
% filenames_new: sorted file names
% MJ
% YW document 10/10/2018
if ~isempty(str2num(filenames{1}(1)))
    %First character in first file is a number
    try
        for i = 1:length(filenames)
            number=strsplit(filenames{i},'_');
            numbers(i)=str2num(number{1});
            [numbers_sorted,idx] = sort(numbers);
            filenames_new=filenames(idx);
        end
    catch exception
        disp('Sort by filename')
        filenames_new = filenames;
    end
else
    filenames_new = filenames;
end
