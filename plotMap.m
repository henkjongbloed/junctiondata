function plotMap(U,S)

%
if U{1}.DS==2015
    latlim = 51.893825 + .006*[-1.5, 1];
    lonlim = 4.319825 + .0025*[-1, 1];
else
    latlim = 51.865964 + .006*[-1, 1];
    lonlim = 4.332862 + .0025*[-1, 1];
end
figure;
geolimits(latlim,lonlim)
geobasemap("satellite")
hold on

for ct = 1:numel(S)
    [lat,lon] = utm2ll(S{ct}.Xold, S{ct}.Yold, 31, 'wgs84');
    geoplot(lat,lon,':w')
end


end