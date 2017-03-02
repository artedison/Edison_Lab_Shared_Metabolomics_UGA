%%
% Written 11_7_13 by Chaevien Clendinen to combine chemical shifts
%AS Edison Lab, University of Florida
%X is a structure 
%fieldname which matrix in the structure X
%range is the ppm flexibility allowed when combining shifts. default is
%0.03

function out=comchem(X,fieldname,range)

if isempty(range)
    range=0.03;
end

for i=1:length(X)
    
    z=zeros(size( X.(fieldname)));
    for j=1:length( X.(fieldname))-1
        if  X.(fieldname)(j)>=X.(fieldname)(j+1)-range &&  X.(fieldname)(j)<= X.(fieldname)(j+1)+range 
            z(j:j+1)=[ X.(fieldname)(j)  X.(fieldname)(j+1)];
        end
    end
    
zi=1;
for j=1:length(X(i).(fieldname))-1
    if z(j)==0
        continue
    elseif z(j)>=z(j+1)-range && z(j)<=z(j+1)+range
        continue
    else
        a=[z(zi:j) 0];
        id=[0 logical(a) 0];
           ii1=strfind(id,[0,1]);
         ii2=strfind(id,[1,0])-1;
        for k=1:length(ii1)
            a(ii1(k):ii2(k))=median(a(ii1(k):ii2(k)));
        end
        z(zi:j)=a(1:end-1);
        zi=j+1;
    end
end

z(zi:end)=median(z(zi:end));


 X.(fieldname)(z~=0)=z(z~=0);
 X.(fieldname)=unique( X.(fieldname));
 out=X.(fieldname);

clear z zi a
end




