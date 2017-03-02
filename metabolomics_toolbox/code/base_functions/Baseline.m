%Baseline Correction for NMR Spectroscopic Metabolomics Data Analysis
%Baseline Model V1.0
%Yuanxin Xi, IDAV, UC Davis
%http://www.biomedcentral.com/1471-2105/9/324
%Input:
%       b Nx1 the spectrum data
%       A 1x1 smoothing factor
%       B 1x1 negative penalty
%       s 1x1 noise standard deviation
%Output:
%       bd Nx1 baseline

function bd=Baseline(b,A,B,s)
    L=length(b);
    Bs=-B/s;  As=-A/s;
    bd=ones(L,1)*median(b);bd0=b;nm=norm(bd-bd0);nm0=realmax;
    M0=-ones(L,1)/As;
    e=ones(L,1);
    D0=spdiags([2*e -8*e 12*e -8*e 2*e],-2:2,L,L);
    D0(1,1)=2; D0(L,L)=2;
    D0(2,1)=-4; D0(1,2)=-4; D0(L,L-1)=-4; D0(L-1,L)=-4;
    D0(2,2)=10;D0(L-1,L-1)=10;
    I=0;
    while nm>10 & I<30
        %& nm<nm0;
        I=I+1;
        M=M0;D=D0;bd0=bd;nm0=nm;
        for i=1:L
            if bd(i)>b(i)
                M(i)=M(i)+2*Bs*b(i)/As;
                D(i,i)=D(i,i)+2*Bs/As;
            end
        end
        bd=D\M;
        nm=norm(bd0-bd);
    end
