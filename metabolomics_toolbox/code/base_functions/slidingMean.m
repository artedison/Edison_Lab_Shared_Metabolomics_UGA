function x = slidingMean(x,window)
%x=[1:10];
%window = 5; %must be odd
startidx = ceil( window / 2 );
endidx = length(x) - floor( window / 2);
range = floor( window / 2);
for i=startidx:endidx
    x(i)=mean(x(i-range:i+range));
end