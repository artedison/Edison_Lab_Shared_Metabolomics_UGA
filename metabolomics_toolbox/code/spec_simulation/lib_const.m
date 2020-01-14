%%% please don't run unless you understand this script.
%%% contact me Yue.Wu@uga.edu if you want to use this script
%%%download gissmo library
%%%load gissmo library into a matlab dataset 600MHZ
%%%spectra between ppm range [-1 11] will be selected
%%%% yue wu
close all;
clear all;
%%wget command
PATH=getenv('PATH');
setenv('PATH',[PATH ':/usr/local/bin/']);
comp='/Users/yuewu/';%the computer user location
%%if redownload is needed do change to another folder
libdir=[comp 'YOURPATH'];%Dropbox (Edison_Lab@UGA)/Resources/gissmo_lib/raw/
%%% downloading
cd(libdir);
tabinfor=readtable([libdir 'infor.txt'],'Delimiter','\t');
%%deal with the lines with simulation_2
wrongind=find(~contains(table2cell(tabinfor(:,5)),'Download'));
tabinfor(wrongind,5)=tabinfor(wrongind,3);
tabinfor(wrongind,3)=tabinfor(wrongind,1);
tabinfor(wrongind,4)=num2cell(str2double(table2cell(tabinfor(wrongind,2))));
tabinfor(wrongind,1:2)=tabinfor(wrongind-1,1:2);
tabinfor(:,2)=strrep(table2cell(tabinfor(:,2)),'BMRB ','');
tabinfor(:,3)=strrep(table2cell(tabinfor(:,3)),'Simu','simu');
compdnames=unique(tabinfor.CompoundName);
%length(compdnames) 1250
homlink='http://gissmo.nmrfam.wisc.edu/entry/';
options=weboptions('Timeout',500);
mkdir(libdir)
for i=1:size(tabinfor,1)
  inforvec=table2cell(tabinfor(i,:));
  weblink=[homlink inforvec{2} '/' inforvec{3} '/zip'];
  diri=num2str(i);
  mkdir(diri);
  cd(diri);
  system(['wget ',weblink]);
  system('unzip zip');
  cd('../');
  pause(5);
end
save('gissmo_tab.mat','tabinfor')

%%load and store the gissmo library
dataloc='/B0s/sim_600MHz';%% the relative location in gissmo folder
peaklistloc='/peaks/sim_600MHz_peaks_standard.csv';
ppmlimit=[-1 11];%the range of ppm that will be inlcuded in the data
% strdataannotbmrb=table2cell(tabinfor(:,2)); %%ids for the compound in order
strdata=struct();
strdatalist=struct();
ppm=[];
%%read data from gissmo folder
for i=1:size(tabinfor,1)
  i
  inforvec=table2cell(tabinfor(i,:));
  %%spectra
  name=[inforvec{2} '_' inforvec{3}];
  % name=strrep(name,'-','_');
  % name=['compd_' strrep(name,'''','_')];
  filelocat=[num2str(i) '/' inforvec{2} '/' inforvec{3} dataloc];
  filelocat=ls([filelocat '*']);%accomodate files with different names
  filelocat=regexprep(filelocat,'[\n\r]+','');
  tab=readtable(filelocat);
  tab=table2array(tab);
  ragind=find(tab(:,1)>ppmlimit(1)&tab(:,1)<ppmlimit(2));
  tab=tab(ragind,:);
  strdata.(name)=tab(:,2);
  ppm=tab(:,1)';
  %%peaklist
  filelocat=[num2str(i) '/' inforvec{2} '/' inforvec{3} peaklistloc];
  tab=readtable(filelocat);
  tab=table2array(tab);
  strdatalist.(name)=tab;%%PPM, amplititude
end
save('wholegissmo_spectral.complex.more.mat','strdata','ppm');
save('wholegissmo_spectral.real.list.mat','strdatalist');
