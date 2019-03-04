function [freq,ps]=philtimeps(time,y)

Fs = 1./(time(2) - time(1)); 
N = length(y);
xdft = fft(y);
xdft = xdft(1:floor(N/2)+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(y):Fs/2;
ps = psdx;

return