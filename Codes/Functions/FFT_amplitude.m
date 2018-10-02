function [y, f] = FFT_amplitude(data, delta)
dt = 1/delta;
N = length(data);
nf = N/2+1;
df = 1/(dt*N);
fmax = 1/dt/2;
f = 0:df:fmax;
ampl_fft = dt*abs(fft(data));
y(1) = ampl_fft(1);
for i=2:nf
    y(i) = ampl_fft(i)+ampl_fft(N-i+2);
end