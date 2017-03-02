function XAL=star_align2D(X,label,ppm1,ppm2,maxshift)

% XAL=star_align2D(X,label,ppm1,ppm2,maxshift)
%
% Star alignment against mean spectrum using local alignment function from HATS.
%
% Arguments:
% 
% X            Data matrix of spectra
% label        label matrix from segment2D.m
% ppm1         Chemical shift vector of F2
% ppm2         Chemical shift vector of F1
% maxshift     Maximum shift if F2 and F1 dimensions [F2,F1]
%
% Return Values:
% XAL          Aligned spectral matrix

XAL=X;

if exist('maxshift')~=1
    if abs(max(ppm1)-max(ppm2))<10
        maxshift=[0.1 0.1];
    else
        maxshift=[1 0.1];
    end 
end
resolution=[size(X,1),size(X,2)];
shift1=round(maxshift(2)/(ppm2(1)-ppm2(2)));
shift2=round(maxshift(1)/(ppm1(1)-ppm1(2)));

%Here we want to avoid aligning a concatenation of spectral peaks, so
%reverse the order
ind=2;
nonind=1;

ref=mean(X,3);
stats=slrregionprops(label);

for spectcount=1:size(X,3)
    wait = waitbar(0,['Aligning ',num2str(spectcount),' of ',num2str(size(X,3))]);
    for reg=1:max(max(label))
        if (stats(reg).Area)/(resolution(1)*resolution(2))<=.001
            region=[floor(stats(reg).BoundingBox(2)),ceil(stats(reg).BoundingBox(2)+stats(reg).BoundingBox(4)), floor(stats(reg).BoundingBox(1)),ceil(stats(reg).BoundingBox(1)+stats(reg).BoundingBox(3))];
            
            
            if region(1)-shift1>=1 && region(2)+shift1<=resolution(1) && region(3)-shift2>=1 && region(4)+shift2<=resolution(2)
                
                for k=-1*shift1:1:shift1 % allow alignment up to align_window threshold ppm
                    for z=-1*shift2:1:shift2
                        reg1=ref(region(1)+k:region(2)+k,region(3)+z:region(4)+z);
                        reg1line=reshape(reg1,[size(reg1,1)*size(reg1,2),size(reg1,3)])';
                        reg2=XAL(region(1):region(2),region(3):region(4),spectcount);
                        reg2line=reshape(reg2,[size(reg2,1)*size(reg2,2),size(reg2,3)])';
                        mat(k+shift1+1,z+shift2+1)=mean(mean(pearson(reg1line,reg2line)));
                    end
                end
                
                if max(max(mat))>=.5
                    [i,j]=ind2sub(size(mat),find(mat==max(max(mat))));
                    %need to take indices from region in moved spectrum (1 in this case) and
                    %move over to optimized corrcoef location (greater than
                    %thresho
                    
                    if (i~=(shift1+1) || j~=(shift2+1))
                        transfer=XAL(region(1):region(2),region(3):region(4),spectcount);
                        
                        %Interpolation at boundaries plus noise
                        if j-(shift2+1)<=-1
                            horz=median(XAL(region(1):region(2),region(4)+(j-(shift2+1)):region(4),spectcount),2);
                        elseif j-(shift2+1)>=1
                            horz=mean(XAL(region(1):region(2),region(3):region(3)+(j-(shift2+1)),spectcount),2);
                        else
                            horz=zeros(region(2)-region(1)+1,1,size(transfer,3));
                        end
                        horzmat=repmat(horz,1,region(4)-region(3)+1);
                        
                        
                        
                        if i-(shift1+1)<=-1
                            vert=mean(XAL(region(2)+(i-(shift1+1)):region(2),region(3):region(4),spectcount),1);
                        elseif i-(shift1+1)>=1
                            vert=mean(XAL(region(1):region(1)+(i-(shift1+1)),region(3):region(4),spectcount),1);
                        else
                            vert=zeros(1,region(4)-region(3)+1,size(transfer,3));
                        end
                        vertmat=repmat(vert,region(2)-region(1)+1,1);
                        regnoise=median(median(XAL(region(1):region(2),region(3):region(4))));
                        r = -.1*regnoise + .2*regnoise.*rand(size(vertmat+horzmat));
                        
                        if max(max(max(abs(horzmat))))>0 && max(max(max(abs(vertmat))))>0
                            for index31=1:size(r,3)
                                prereplace{index31}=cat(3,horzmat(:,:,index31),vertmat(:,:,index31));
                                [c,location]=min(abs(prereplace{index31}),[],3);
                                for d1=1:size(location,1)
                                    for d2=1:size(location,2)
                                        replace2(d1,d2,index31)=prereplace{index31}(d1,d2,location(d1,d2));
                                    end
                                end
                                clear c location
                            end
                            replace=replace2+r;
                        elseif max(max(max(abs(horzmat))))>0 && max(max(max(abs(vertmat))))==0
                            replace=horzmat+r;
                        else
                            replace=vertmat+r;
                        end
                        
                        XAL(region(1):region(2),region(3):region(4),spectcount)=replace;%replace region with interpolation
                        XAL(region(1)+i-(shift1+1):region(2)+i-(shift1+1),region(3)+j-(shift2+1):region(4)+j-(shift2+1),spectcount)=transfer; %this shifts the real region
                        clear replace replace2 transfer transfer2 horz horzmat vert vertmat r regnoise
                    end
                end
            end
             waitbar(reg / max(max(label)))
        end
    end
    close(wait)
    
end
return;
end


function stats=slrregionprops(label)

for k=1:max(max(label))
    
    [dim(:,1),dim(:,2)] = find(label==k);
    stats(k).Area=size(dim,1);
    stats(k).Centroid=mean(dim);
    start=min(dim,[],1);
    boxlength=(max(dim,[],1)-min(dim,[],1))+1;
    stats(k).BoundingBox(1)=start(2)-.5;
    stats(k).BoundingBox(2)=start(1)-.5;
    stats(k).BoundingBox(3)=boxlength(2);
    stats(k).BoundingBox(4)=boxlength(1);
    clear dim
    
end

stats=stats';

return

end


function [R, P, Rlo, Rup] = pearson(A,B)

% PEARSON computes the Pearson (linear) correlation matrix.
% ----------------------------
% [R P Rlo Rup] = pearson(A,B)
% ----------------------------
% Description: Computes the Pearson (linear) correlation matrix. There are
%              two differences between this function and CORRCOEF. First,
%              the data is variables-by-samples here, in contrast to
%              samples-by-variables in CORRCOEF. Second, corrcoef(A,B) is
%              the same as corrcoef([A B]). This is not the case for
%              PEARSON, where pearson([A ; B]) = corrcoef([A' B']) but
%              pearson(A,B) is just a (potentially rectangular) matrix of
%              the correlation for each pair (a,b) where a is in A and b is
%              in B.
% Input:       {A} any data matrix variables-by-samples.
%              <{B}> any data matrix with the same number of
%                   measurement as in {A}. Basically, the correlation
%                   coefficients are calculated between the variables of
%                   {A} and {B}. If {B} is missing, the correlation
%                   coefficients are calculated between the variables of
%                   {A}.
% Output:      {R} Pearson correlations.
%              {P} P-value of the independence hypothesis testing.
%              {Rlo} lower bound for 95% confidence interval.
%              {Rup} upper bound for 95% confidence interval.

% ? Liran Carmel
% Classification: Correlations
% Last revision date: 20-Mar-2006

% low-level input parsing
error(nargchk(1,2,nargin));
single_matrix_correlation = false;
if nargin == 1
    single_matrix_correlation = true;
end

% calculate
if single_matrix_correlation
    [R P Rlo Rup] = corrcoef(A');
elseif isscalar(A)
    R = nan;
    P = 1;
    Rlo = nan;
    Rup = nan;
else
    vars_A = 1:size(A,1);
    vars_B = size(A,1) + (1:size(B,1));
    [R P Rlo Rup] = corrcoef([A' B']);
    R = R(vars_A,vars_B);
    P = P(vars_A,vars_B);
    Rlo = Rlo(vars_A,vars_B);
    Rup = Rup(vars_A,vars_B);
end
end

