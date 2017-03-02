%############################################################################
%#
%#                          function PEAKPICK
%#
%#      finds local maxima (peaks)
%#
%#      usage: peaks=peakpick(data,fw,offset,plot);
%#                 
%#                   data = input spectrum (1D, real)
%#                  peaks = output index of peak found
%#                     fw = OPTIONAL, smoothing strength, default=20
%#                 offset = OPTIONAL, noise offset (green line in plot)
%#                   plot = OPTIONAL, if not 0 plot is generated
%#              
%#              example: pp=peakpick(data,500,0.01*max(data),1);
%#                       data(pp)     
%#         
%#      (c) P. Blümler 1/03
%############################################################################
%----------------------------------------------------------------------------
%  version 1.1 PB 21/1/03    (please change this when code is altered)
%----------------------------------------------------------------------------
function peaks=peakpick(data,fw,offset,plot_flag)

data=squeeze(data);
% dim=dimension(data);
% if dim ~= 1
%     errordlg('ERROR: data input array is NOT ONEDIMENSIONAL!');
%     return
% end

switch nargin
case 1, fw=20; 
        offset=3*std(data(round(0.9*end):end));
        plot_flag=0;
case 2, offset=3*std(data(0.9*end:end));
        plot_flag=0;  
case 3, plot_flag=0;  
end

%smooth data
x=linspace(-1,1,length(data));
filter=exp(-x.^2*fw);
temp_data=fftshift(fft(data)).*filter;
temp_data=real(ifft(ifftshift(temp_data)));

dd1=gradient(temp_data);
dd2=gradient(dd1);
sdd1=abs(gradient(sign(dd1)));
sdd2=abs(gradient(sign(dd2)));
peaks=find(sdd1 ~=0 &  dd2 < 0 & sdd2==0 & data>offset);
for t=2:length(peaks)
  if peaks(t)==peaks(t-1)+1
       peaks(t-1)=-1;
    end
end
peaks=sort(peaks);
index=min(find(peaks~=-1));
peaks(1:index-1)=[];

if plot_flag~=0
  plot(data);hold on;
  plot(peaks,data(peaks),'ro');
  plot([1:length(data)],offset,'g','Linestyle','-.');
  hold off
end
