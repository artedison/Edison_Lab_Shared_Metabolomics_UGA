function [TP,F1,FP,FN] = getF1(compounds, peaks)
%% example
% peaks = [14.*rand(100,1), 140.*rand(100,1)];
% 
% compoundnames = {'Leucine','meucine','neucine'};
% compounds{1} = [peaks([4,34,23],:);peaks(5,:)+[.02,-.02];[34,210]];
% compounds{2} = [peaks([5,36,23,67],:);[10,210];[30,210]];
% compounds{3} = [peaks([6,35],:);[10,200];[24,200]];
%%
for i=1:length(compounds)
[matches,TP(i),FP(i),FN(i)] = intersectTol(compounds{i},peaks,[.05,.15]);

Precision = TP(i)/(TP(i)+FP(i));
Recall = TP(i)/(TP(i)+FN(i));
F1(i) = 2*(Precision*Recall)/(Precision+Recall);
end

end

function [matches,TP,FP,FN] = intersectTol(A, B, tol)

matches = [];
FP=0;
idxA=zeros(size(A,1),1);
idxB=zeros(size(B,1),1);
for i=1:size(A,1)
    Cpeak = A(i,2);
    idx = abs(Cpeak-B(:,2))<tol(2);
    Bmatches = B(idx,:);
    for j=1:size(Bmatches,1)
        FP = FP + size(Bmatches,1)-1;
        if all(abs(A(i,:)-Bmatches(j,:))-tol<=0)
            idxA(i)=1;
            idxB(i)=1;
            matches = [matches;[A(i,:), Bmatches(j,:)]];
        end
    end
end
TP=sum(idxA);
FN=sum(~idxA);
end