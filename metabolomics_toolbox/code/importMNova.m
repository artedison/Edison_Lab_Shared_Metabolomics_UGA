function importMNova(filename)
if isa(filename,'char')
    A=importdata(filename); %one file containing many 1ds or one 2d or one 1d
else
    A=importdata(filename{1}); %multiple files each containing one 1d or one 2d
end
norm=1;
if isa(filename,'char') & isnumeric(A) & size(A,2)==2 %one file containing one 1d
    X=A(:,2)';
    ppm=A(:,1)';
    assignin('base','X',X)
    assignin('base','ppm',ppm)
    assignin('base','XTitles',filename)
elseif isa(filename,'char') & isnumeric(A) & size(A,2)~=2 %one file containing one 2d
    spectra.ppm1=A(1,[2:end]);
    spectra.ppm2=A([2:end],1);
    spectra.real=A([2:end],[2:end]);
    assignin('base','XTitles',filename)
    assignin('base','spectra',spectra)
    norm=0;
elseif isa(filename,'char') & ~isnumeric(A) %one file containing multiple 1ds
    X=A.data(2:end,:)';
    ppm=A.data(1,:);
    assignin('base','X',X)
    assignin('base','ppm',ppm)
    assignin('base','XTitles',A.textdata(2:end,2:end))
elseif isa(filename,'cell') & size(A,2)==2 %multiple files each containing one 1d
    for i=1:length(filename)
        A=importdata(filename{i});
        spectra(i).real=A(:,2)';
        spectra(i).ppm=A(:,1)';
        spectra(i).Title=filename{i};
    end
    [X,ppm,XTitles]=Setup1D(spectra);
    assignin('base','X',X)
    assignin('base','ppm',ppm)
    assignin('base','XTitles',XTitles)
elseif isa(filename,'cell') & size(A,2)~=2 %multiple files each containing one 2d
    for i=1:length(filename)
        A=importdata(filename{i});
        spectra(i).ppm1=A(1,[2:end]);
        spectra(i).ppm2=A([2:end],1);
        spectra(i).real=A([2:end],[2:end]);
        spectra(i).Title=filename{i};
    end
    Setup2D(spectra)
    assignin('base','XTitles',XTitles)
    assignin('base','X',X)
    assignin('base','XNoise',XNoise)
    assignin('base','ppm1',ppm1)
    assignin('base','ppm2',ppm2)
    norm=0;
else
    error('Unrecognized file, please check format')
end

if norm
    %normalize by height of TSP
    [~,idx]=min(abs(ppm));
    tsp=X(:,idx);
    nf=repmat(tsp/max(tsp),1,size(X,2));
    assignin('base','X_tsp',X./nf)
end

%normalize by area of TSP
% [~,idx1]=min(abs(ppm-.005));
% [~,idx2]=min(abs(ppm+.005));
% tsp2=sum(X(:,idx2:idx1),2);
% nf2=tsp2/max(tsp2);
% X_tsp2=X./nf2;

%fft interpolation
% ppm2=linspace(min(ppm),max(ppm),32768);
% xi=interpft(x,32768);


