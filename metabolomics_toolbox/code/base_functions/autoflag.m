function flag_descriptor=autoflag(Y,Yvec,XTitles,Xtitles,filename)
%%
% Y should be a vector of the parameter you want tested. 
% Yvec identity vector defining the groups
% XTitles are the samople id in vector form
% Xtitles are the sample ids in cell array
% filename should be in quotes to make it a string. This is what the ame of
%           the figure filename will be. 
%This will standardize the data using y-mu/stan dev.
% this will then create a new vector named as flag_descriptor which will
% contain zeros for within 1 sd and flagged as a 1 if outsie of 1sdevs.
% Can change to be however many stdevs you want, but NMR data has shown to
% be pretty reproducible so more stirct seems to be ok. 


% Disclaimer: the Yvec should contain no zeros. 

% Brittany Lee
% 2015, UF


id=[XTitles Y];
groups=unique(Yvec);
flag_descriptor=[XTitles zeros(size(Y))];

for k=1:size(Y);
    for i=1:length(groups);
        for j=groups(i);
            if Yvec(k)==j;
               y=Y(k);
               M=mean(Y(Yvec==i)) ;% if you get 'Subscript indices must either be real positive integers or logicals.'
               SD=std(Y(Yvec==i)); % then need to change the zeros in Yvec to be non zero numbers
               z=((Y(k)-M)/SD);
              if abs(z)>1==1;  % change here how big or small youw ant the z to be within
                flag_descriptor(k,2)=1;
               else
                   break
              end
             else Yvec(k)~=j;
               break
            end
            end
    end
    end
       
h=figure
hold on
scatter(Yvec(flag_descriptor(:,2)==1),[Y(flag_descriptor(:,2)==1)],'r')
scatter(Yvec(flag_descriptor(:,2)==0),[Y(flag_descriptor(:,2)==0)],'b')
text(Yvec+0.02,Y,Xtitles,'FontSize', 15)
title('Distribution per group - Flagged in Red','FontSize', 18)
    xlabel('Group')
    ylabel('AU')
hold off


print (h, '-dpdf',filename)
print (h, '-depsc',filename)

