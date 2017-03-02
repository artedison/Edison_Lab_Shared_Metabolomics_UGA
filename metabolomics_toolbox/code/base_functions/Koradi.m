function [noise,basenoise]=Koradi(a)

warning off
for z1=1:size(a,1)
    for index=1:16
        stdmatdir(z1,index)=std(a(z1,((index-1)/16)*size(a,2)+1:(((index/16)*size(a,2))-1)));
    end
end
stddir=min(stdmatdir,[],2);

for z1=1:size(a,2)
    for index=1:16
        stdmatindir(index,z1)=std(a(((index-1)/16)*size(a,1)+1:(((index/16)*size(a,1))-1),z1));
    end
end
stdindir=min(stdmatindir,[],1);

basenoise=min([min(stdindir,[],2), min(stddir,[],1)]);
addstdindir=sqrt((stdindir.^2)-basenoise^2);
addstddir=sqrt((stddir.^2)-basenoise^2);

noise=zeros(size(a,1),size(a,2));
for z=1:size(a,1)
    noise(z,:)=sqrt((addstddir(z,:)^2)+(addstdindir.^2)+basenoise^2);
end