function Topnames=Topnames(Plist,fieldname)

for q=1:length(Plist)
    for j=1:length(Plist(q).(fieldname))
        for i=1:length(Plist)
            names(i,:)= strcmp(Plist(q).(fieldname)(j),Plist(i).(fieldname));
            if sum(sum((fieldname)))>=3
                Topnames(q,j)=Plist(q).(fieldname)(j);
            end
        end
        clear names
    end
end

Topnames=reshape(Topnames(~cellfun(@isempty,Topnames)),1,[])';

Topnames=unique(strcat(Topnames(:,1)));
end