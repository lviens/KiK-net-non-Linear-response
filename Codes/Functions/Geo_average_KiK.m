function [aver, std1, std2] = Geo_average_KiK(data)

% Compute geometric mean
aver = geomean(data);

for k = 1:length(data)
    si = 0;
    for l = 1:length(data(:,1))
        si = (log(data(l,k)) - log(aver(k)))^2 + si;
    end    
    % One std to the mean
    std1(1,k) = exp(log(aver(k)) + sqrt(si/(length(data(:,1))-1)));
    std1(2,k) = exp(log(aver(k)) - sqrt(si/(length(data(:,1))-1)));
    
    %  Two std to the mean
    std2(1,k) = exp(log(aver(k)) + 2*sqrt(si/(length(data(:,1))-1)));
    std2(2,k) = exp(log(aver(k)) - 2*sqrt(si/(length(data(:,1))-1)));
end