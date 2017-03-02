function write_to_txt_v2(filename, header, left_columns, data, sep)
% write_to_txt(filename, header, left_columns, data, sep)

if ~exist('sep')
    sep = char(9);
end

if ~isempty(left_columns) && ~isempty(data)  % check if data and left_columns are consistent
    if size(data,1)~=size(left_columns,1)
        error('row numbers do now match between left_columns and data');
    end
end


fid = fopen(filename,'w');
% define and write header
if ~isempty(header)
    header_row = cell2mat(strcat(header,{sep})); header_row(end)=char(10);
    fwrite(fid,header_row,'char'); 
end

if isempty(left_columns) && isempty(data) % if only header exists
    fclose(fid);
    return
end

total_columns = size(left_columns,2)+size(data,2);
total_rows = max(size(left_columns,1),size(data,1));

block_size = 1000000;
row_block_size = ceil(block_size/total_columns);
for i=1:row_block_size:total_rows
    % cover data into left_column-like
    if isempty(data)
        data_columns=[];
    else
        data_tmp = data(i:min(size(data,1),i+row_block_size-1),:);
        num_columns = size(data_tmp,2);
        num_rows = size(data_tmp,1);
        data_tmp_str = num2str(reshape(data_tmp,1,num_rows*num_columns));
        data_tmp_str(regexp(num2str(data_tmp_str),' +','start'))=sep;
        data_tmp_str(data_tmp_str==' ')=[];
        data_tmp = reshape(regexp(data_tmp_str,char(9),'split'),num_rows,num_columns);
        clear data_tmp_str
    end
    
    if isempty(left_columns)
        left_columns_tmp=[];
    else
        left_columns_tmp = left_columns(i:min(size(data,1),i+row_block_size-1),:);
    end

    % merge left_column and data_column and write
    this_block_to_write = [left_columns_tmp, data_tmp];
    num_columns = size(this_block_to_write,2);
    num_rows = size(this_block_to_write,1);
    this_block_to_write_str = cell2mat(reshape(strcat(this_block_to_write,{sep})',1,num_columns*num_rows));
    ind = find(this_block_to_write_str==sep);
    this_block_to_write_str(ind(num_columns:num_columns:end)) = char(10);
    fwrite(fid,this_block_to_write_str,'char');
    clear this_block_to_write this_block_to_write_str left_columns_tmp data_tmp
end
fclose(fid);
return