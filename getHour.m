function Hour = getHour(t, Trans)

Hour = NaN(size(t));
hourloc = 1;
Hour(1) = hourloc;
for i = 2:length(t)
    if t(i)-t(i-1) > 4500
        hourloc = hourloc + 1;
    end
    if Trans(i)-Trans(i-1) > 0
        hourloc = 1;
    end
    Hour(i) = hourloc;
end
        
