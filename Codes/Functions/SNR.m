
function [r] = SNR(signal,noise)
%  Compute signal-to-noise ratio 
aps = mean(signal.^2); % average power of signal
apn = mean(noise.^2);  % average power of noise
r = 10*log10(aps/apn); % signal-to-noise ratio in decibel
return