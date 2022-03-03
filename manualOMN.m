function t = manualOMN(told, adcpH)
% Purpose: Manually change signal told such that it fits adcpH in terms of
% values.

t2 = rescale(told, min(adcpH), max(adcpH));
hour = datenum('12-Aug-2014 01:00:00') - datenum('12-Aug-2014 00:00:00');
% plot(adcpH)
% hold on
% plot(t2)
uh = unique(adcpH);
t3 = NaN(size(told));
for i = 1:length(uh)
    ind = abs(t2 - uh(i)) < hour/2;
    t3(ind) = uh(i);
end


t = t3;
% plot(t3)


end