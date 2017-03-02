function Total=Loop_Threshes(compounds,S)

CovarThreshold=[0.3 0.4 0.5 0.6 0.7 0.8 0.9];
PvalThreshold=[0.01 0.02 0.03 0.04 0.05];
options=[1 2];
FDR=[0 1];
ShiftRange=[0.03,0.03;0.05,0.05;0.08 0.08;0.1 0.1];

tol=[0.5,0.10];
k=0;
for a=1:length(CovarThreshold)
    for b=1:length(PvalThreshold)
        for c=1:length(options)
            for d=1:length(FDR)
                for e=1:length(ShiftRange)
                    [~,PairsX2H]=genPeaklist(S.corrHXpp,S.covarHXpp,10,S.ppmHpp,S.ppmXpp,'CovarThreshold',CovarThreshold(a),'PvalThreshold',PvalThreshold(b),'options',options(c),'FDR',FDR(d),'ShiftRange',ShiftRange(e,:));
                    GrandList=1;
                    h=1;
                    for j=1:length(PairsX2H)
                        if isempty(PairsX2H(j).toppm)
                            continue
                        end
                        GrandList(h:h+length(PairsX2H(j).toppm)-1,1)=PairsX2H(j).toppm;
                        GrandList(h:h+length(PairsX2H(j).toppm)-1,2)=PairsX2H(j).ppm_of_carbon;
                        h=length(GrandList)+1;
                    end
                    
                    for n=1:length(compounds)
                        A=compounds{n};
                        
                        TPc=[];
                        FPc=[];
                        FNc=[];
                        TNc=[];
                        
                        for f=1:size(A,1)
                            Cpeak = A(f,2);
                            idx = abs(Cpeak-GrandList(:,2))<tol(2);
                            Bmatches = GrandList(idx,:);
                            zz=length(find(A==0));
                            
                            TPtemp=0;
                            FPtemp=0;
                            for j=1:size(Bmatches,1)
                                if  sum(abs(Bmatches(j,1)-A(:,1))<=tol(1))~=0
                                    TPtemp=TPtemp+1;
                                else
                                    FPtemp=FPtemp+1;
                                end
                            end
                            TPc(f)=TPtemp;
                            FPc(f)=FPtemp;
                            FNc(f)=size(A,1)-zz-TPc(f);
                            TNc(f)=length(GrandList)-size(A,1)-zz;
                        end
                        
                        TPi(n)=mean(TPc);
                        FPi(n)=mean(FPc);
                        FNi(n)=mean(FNc);
                        TNi(n)=mean(TNc);
                    end
                    k=k+1;
                    Total(k).set=[CovarThreshold(a) PvalThreshold(b) options(c) FDR(d) ShiftRange(e,r)];
                    Total(k).TP=sum(TPi);
                    Total(k).FP=sum(FPi);
                    Total(k).FN=sum(FNi);
                    Total(k).TN=sum(TNi);
                    
                    Total(k).Sensitivity=Total(k).TP/(Total(k).TP+Total(k).FN);
                    Total(k).Specificity=Total(k).TN/(Total(k).FP+Total(k).TN);
                end
            end
        end
    end
end

end
