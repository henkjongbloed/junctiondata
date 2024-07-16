function [T, ind] = getHour(t, dt)
% Detects sudden changes in time series (essentially a moving average
% filter with variable window size...

T = NaN(size(t));
ind = NaN(size(t));
% hr_begin = 1;
idx_begin = 1;
idx = 0;
for i = 1:length(t)-1
    if (t(i+1)-t(i) > dt) || (i == length(t) - 1)
        idx_end = i;
        T(idx_begin:idx_end) = mean(t(idx_begin:idx_end));
        ind(idx_begin:idx_end) = idx;
        idx = idx + 1;
        idx_begin = idx_end + 1;
    end
end

T(end) = T(end-1);
ind(end) = ind(end-1);
% figure;
% plot(T)
% title('getHour')
% subplot(211)
% plot(t)
% hold on
% plot(T)
% subplot(212)
% plot(ind)

end