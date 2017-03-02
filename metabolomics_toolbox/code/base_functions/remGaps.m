function out=remGaps(input)
%remove gaps in the numbering
%out=remGaps([0,0,0,1,1,8]) -> [1,1,1,2,2,3]

a=1; 
b=1;
df2=zeros(size(input));
dfl=length(input);
while a<=dfl
    df2(input==input(a))=b;
    input(input==input(a))=0;
    b=b+1;
    a=find(input,1,'first');
end
out=df2;