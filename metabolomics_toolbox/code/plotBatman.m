function plotBatman()
% plot Batman fitted spectra from specFit_*.txt
% GSS

txtlist = dir('specFit_*.txt');
txtlist = {txtlist.name};
txtlist=cellfun(@(x) ([cd filesep x]),txtlist,'uniformoutput',0);

for i=1:length(txtlist)
  BX=importSpecFit(txtlist{i});
  ppm=BX(:,1);
  original_spectrum = BX(:,2);
  metabolite_fit = BX(:,3);
  wavelet_fit = BX(:,4);
  overall_fit = BX(:,5);
  
  figure;
  h1=subplot(2,1,1);
  hold on;
  plotr(ppm,original_spectrum,'b');
  plotr(ppm,metabolite_fit,'r');
  legend('Original Spectrum','Metabolite Fit')
  hold off;
  h2=subplot(2,1,2);
  hold on;
  plotr(ppm,original_spectrum,'b');
  plotr(ppm,wavelet_fit,'g');
  legend('Original Spectrum','Wavelet Fit')
  %plotr(ppm,overall_fit,'k');
  hold off;
  linkaxes([h1,h2],'xy');
  pause;
end
        
        
        
        
