function subplot_with_k_panels(i,k)
if k<=3
    subplot(1,k,i);
elseif k<=8
    subplot(2,ceil(k/2),i);
elseif k<=15
    subplot(3,ceil(k/3),i);
else
    subplot(4,ceil(k/4),i);
end
