function plotMap(mesh, DS)

if DS==2015
    lat0 = 51.8938250;
    lon0 = 4.3198250;
    
else
    lat0 = 51.8659640;
    lon0 = 4.332862;
end
[xr, xl] = utm2ll(mesh.xw, mesh.yw, 31, 'wgs84');
% figure;
% subplot(1,2,1)
makefigure(10,7.5)

latlim = lat0 + 6*.006*[-.5, 4];
lonlim = lon0 + 10*.0025*[-15, 1];
geolimits(latlim ,lonlim)
geobasemap("satellite")
hold on
geoplot(xr, xl,'w*', 'MarkerSize',5)
% struct(gca)
% set(gca, 'LongitudeString', 'Longitude', 'interpreter', 'latex')
% set(gca, 'LatitudeString', 'Latitude', 'interpreter', 'latex')

% gca.XLabel.Interpreter = 'latex';
% gca.YLabel.Interpreter = 'latex';
makefigure(10,7.5)


latlim = lat0 + .006*[-1, 1];
lonlim = lon0 + + .0025*[-1, 1];
geolimits(latlim , lonlim )
geobasemap("satellite")
hold on
geoplot(xr, xl,'w')
geoplot(xr, xl,'w*', 'MarkerSize',5)
XLabel.Interpreter = 'latex';
YLabel.Interpreter = 'latex';
% for ct = 1:numel(S)
%     [lat,lon] = utm2ll(S{ct}.Xold, S{ct}.Yold, 31, 'wgs84');
%     geoplot(lat,lon,':w')
% end


end