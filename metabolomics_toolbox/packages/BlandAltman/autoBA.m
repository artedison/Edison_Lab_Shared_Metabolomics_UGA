function autoBA(Xrand,XTitles) 

%% 
% X should be a peak picked matrix for a particular set of data (i.e all
% pooled samples from an experiment)


%%

for i=1:(size(Xrand,1)-1)
    for b=1:(size(Xrand,1)-1)
        a=i+b;
        if a<=size(Xrand,1)
            labels={XTitles(i),XTitles(a)};
            BlandAltman(Xrand(i,:)',Xrand(a,:)',labels);
        else
           break
        end
        clear labels
    end
    clear b a 
end

