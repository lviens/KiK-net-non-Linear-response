function y = BP_butter_Filter(data, order, delta, cut, cut2)

data = detrend(data);
data = data-mean(data);

w=[cut/(delta/2) cut2/(delta/2)];
[z, p, k]=butter(order,w,'bandpass');
[sos,g]=zp2sos(z,p,k);

y = filtfilt(sos,g,data);
