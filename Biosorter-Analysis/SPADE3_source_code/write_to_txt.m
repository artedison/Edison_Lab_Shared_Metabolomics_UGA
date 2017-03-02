function write_to_txt(filename, header, left_columns, data, sep)
% write_to_txt(filename, header, left_columns, data, sep)

if ~exist('sep')
    sep = char(9);
end


fid = fopen(filename,'w');
% 1st row
if ~isempty(header)
    tmp = header{1};
    for i=2:length(header)
        tmp = [tmp, sep, header{i}];
    end
%     tmp = [tmp, char(13),char(10)];
    fprintf(fid,'%s\n',tmp);
end

for i=1:max(size(left_columns,1),size(data,1))
    i
    tmp = [];
    if ~isempty(left_columns) 
        if isnumeric(left_columns{i,1})
            left_columns{i,1} = num2str(left_columns{i,1});
        end
        tmp = left_columns{i,1};
        for j=2:size(left_columns,2)
            if isnumeric(left_columns{i,j})
                left_columns{i,j} = num2str(left_columns{i,j});
            end
            tmp = [tmp, sep, left_columns{i,j}];
        end
    end
    for j=1:size(data,2)
        if isempty(tmp)
            tmp = [tmp, num2str(data(i,j))];
        else
            tmp = [tmp, sep, num2str(data(i,j))];
        end
    end
    fprintf(fid,'%s\n',[tmp,char(13)]);
end
fclose(fid);
return