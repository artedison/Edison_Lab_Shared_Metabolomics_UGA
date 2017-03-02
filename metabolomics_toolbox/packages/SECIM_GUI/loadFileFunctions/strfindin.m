function idx=strfindin(name,substr)
if iscell(substr)
    idx=find(~cellfun(@isempty,strfind(substr,name)));
else
    idx=find(~cellfun(@isempty,strfind(name,substr)));
end