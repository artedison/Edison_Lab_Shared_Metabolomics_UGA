%% this program calculate simualte surface under different parameter A, B, and C
%% under a power 2 polynomial
%% the H and K curvature will be calculted afterwrads for the middle point(x,y)=(0,0)
%% (x,y,z) z=Ax^2+Bxy+Cy^2
%% expectation: H=A+C K=4AC-B^2
%% as the data is produced from a power 2 polynomial, the orthongonal polynomial should calcualte H and K exactly as expected if the scale/unit is the same
path='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/result/test.surf.plot/';
para=[1 1 -4;%H=-3 K=-17 saddle ridge
    1 1 -1;%H=0 K=-5 minimal surface
    2 6 2;%H=4 K=-20 saddle valley
    -1 2 -1;% H=-2 K=0 ridge surface
    0 0 0;% H=0 K=0 flat surface
    3 6 3;% H=6 K=0 valley surface
    -1 1 -1;% H=-2 K=3 peak surface
    1 1 1;% H=2 K=3 pit surface
    ];%%%A,B,C no condition on K>0 and H=0
namesvec={'Saddle Ridge' 'Minimal Surface' 'Saddle Valley' 'Ridge' 'Flat' 'Valley' 'Peak' 'Pit'};
% [X,Y]=meshgrid(-1:0.02:1,-1:0.02:1);
[X,Y]=meshgrid(-50:1:50,-50:1:50);
windsize=7;
halfwidplus=(windsize+1)/2;
halfwidminus=(windsize-1)/2;
ind=51;%%where the H and K is calculated
for i=1:size(para,1)
  Apara=para(i,1);
  Bpara=para(i,2);
  Cpara=para(i,3);
  Z = Apara*X.^2 + Bpara*X.*Y + Cpara*Y.^2;
  %check whether the measured H and K is the same as expected
  x=1:windsize;%/windsize
  x=x-halfwidminus-1;
  miu=[];
  miu(1)=windsize;%miu0
  miu(2)=0;
  miu(3)=sum(x.^2);%miu2
  phi=[];
  phi(1,:)=repmat(1,1,windsize);
  phi(2,:)=x;
  phi(3,:)=x.^2-repmat(miu(3)/miu(1),1,windsize);
  b=[];
  for j=1:3
    % b(j,:)=phi(j,:)./sqrt(sum(phi(j,:).^2));
    b(j,:)=phi(j,:)./sum(phi(j,:).^2);
  end
  % for j=1:length(x)
  %   b(:,j)=phi(:,j)./sqrt(sum(phi(:,j).^2));
  %    % b(j,:)=phi(j,:)./sum(phi(j,:).^2));
  % end
  indrange=(ind-halfwidminus):(ind+halfwidminus);
  tempmat=Z(indrange,indrange);
  a=zeros(3,3);
  for pi=1:3
    for pj=1:3
      summat=(b(pi,:)'*b(pj,:))'.*tempmat;
      a(pi,pj)=sum(summat(:));
    end
  end
  %f fu fv fuu fvv fuv fvu
  % yfuncvec=[tempmat(halfwidplus,halfwidplus) a(2,1) a(1,2) a(3,1) a(1,3) a(2,2) a(2,2)];
  % a=a';
  yfuncvec=[tempmat(halfwidplus,halfwidplus) a(2,1) a(1,2) 2*a(3,1) 2*a(1,3) a(2,2) a(2,2)];
  % yfuncvec=[0 0 0 2*Apara 2*Cpara Bpara Bpara];
  %% vec: x xu xv xuu xvv xuv xvu the first row not used.
  surfvecmat=[0 0 yfuncvec(1); 1 0 yfuncvec(2); 0 1 yfuncvec(3); 0 0 yfuncvec(4); 0 0 yfuncvec(5); 0 0 yfuncvec(6); 0 0 yfuncvec(7)];%% the first row is never used in calculation
  nv=cross(surfvecmat(2,:),surfvecmat(3,:));
  nval=nv/sqrt(sum(nv.^2));
  G=[surfvecmat(2,:)*surfvecmat(2,:)'  surfvecmat(2,:)*surfvecmat(3,:)'; surfvecmat(2,:)*surfvecmat(3,:)'  surfvecmat(3,:)*surfvecmat(3,:)'];
  B=[surfvecmat(4,:)*nval' surfvecmat(6,:)*nval'; surfvecmat(6,:)*nval' surfvecmat(5,:)*nval'];
  detG=det(G);
  invG=[G(2,2) -G(1,2); -G(2,1) G(1,1)]/detG;
  H=trace(invG*B)/2;
  K=det(B)/detG;
  fig=figure(), hold on
    surf(X,Y,Z,'FaceColor','Interp');
    shading interp;
    set(gca,'box','off');
    title([namesvec{i} ' calcualted: H=' num2str(H) ';K=' num2str(K)]);
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'ztick',[]);
    set(gca,'zticklabel',[]);

  % fig=gcf;
  saveas(fig,strcat(path,'surf.A',num2str(Apara),'.B',num2str(Bpara),'.C',num2str(Cpara),'.fig'));
  close(fig);
end

% for i=3:5
%   for j=3:5
%     cc=a.*(phi(:,i)*phi(:,j)');
%     sum(cc(:))
%   end
% end
cd(path);
filetab=dir('*.fig');
for i=1:size(filetab,1)
  figname=filetab(i).name;
  inputfile=[path figname];
  outputfile=[path strrep(figname, '.fig', '.pdf')];
  fig=openfig(inputfile);
  print(fig,[outputfile],'-dpdf','-fillpage','-r2400');%'-painters','-fillpage',
  close(fig);
end
