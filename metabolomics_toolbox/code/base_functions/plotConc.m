function plotConc()
%plot relative concentration from Batman Output
%
Conc=importdata('RelCon.txt');
XC=Conc.data;
for i=1:size(XC,1)
 figure;
 bar(XC(i,:));
 title(Conc.textdata{i+1,1});
 pause;
end
