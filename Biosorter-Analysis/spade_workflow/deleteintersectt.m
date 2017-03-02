function [data] = deleteintersectt(data,intersectt)
%intersectt = intersectt';
for z=1:length(intersectt);
    data(intersectt{z,1}) = []; 
end