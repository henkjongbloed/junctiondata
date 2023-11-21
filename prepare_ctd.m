function S = prepare_ctd(U, ctd, DS)

S.S          =  ctd.S; %salinity
hour = datenum('12-Aug-2014 01:00:00') - datenum('12-Aug-2014 00:00:00');
% WGS2UTM ook bij ADCPData gebruikt?? ->
if strcmp(DS,'NMOMNW15')    % Change ADCP coordinates rather than CTD coordinates!!
    [S.Xold, S.Yold]      = wgs2utm(ctd.Y, ctd.X); %Convert LatLon to UTM - pay close attention to coordinate systems. WGS84 - ETRS89 differs at most 65 cm, not relevant
    % Time
    [ctdH, S.i] = getHour(ctd.T2, hour/2);% + datenum('9-Aug-2014 00:00:00'); % Correct
    S.t = rescale(ctdH, min(U.adcpH), max(U.adcpH)); % Everything works perfectly
    % Space
    xnew = center_cloud([S.Xold, S.Yold], U.xs.origin'); % Because of offset in CTD coordinates relative to ADCP - quick fix
    S.X = xnew(:,1); S.Y = xnew(:,2); % Only needed for 2015 data

elseif strcmp(DS,'OMHA14')
    [templon, templat, ~]      = rd2wgs(ctd.X,ctd.Y); %Convert RD to LatLon %Correct
    [S.X, S.Y]      = wgs2utm(templat, templon); %Convert LatLon to UTM %Correct
    S.Xold = S.X;
    S.Yold = S.Y;
    ctddt = 6000;
    [ctdH, S.i] = getHour(ctd.T, ctddt); % dt = 6000 -> 1 sample per 6s?
    S.t = rescale(ctdH, min(U.adcpH), max(U.adcpH));
    if strcmp(U.BN, ' Old Meuse North')
        [S.t, S.i] = manualOMN(ctd.T, U.adcpH);
    end

%      max(adcpT) - min(adcpT)
%      max(ctd.T) - min(ctd.T)
%     % Continue here: Make S.t corresponding to 
%     S.t = ctd.T;% + datenum('12-Aug-2014 00:00:00');
    %U.mesh_mean.plot()
end

S.eta = U.eta(S.t); % Waterlevel per CTD sample

S.Z          =  - ctd.Z + S.eta;

[S.s, S.n] = U.xs.xy2sn(S.X, S.Y);

S.zb = n2zb(S.n, U.mesh_mean);
S.sig = z2sig(S.Z, S.zb, S.eta); % IS THIS ALLOWED IN SIGMA COORDINATES?? - TAKE CARE
S.Z = sig2z(S.sig, S.zb, S.eta);

[S.Xp, S.Yp] = U.xs.sn2xy(0*S.s, S.n); % Project on cross section by setting s = 0. [Xp, Yp] = projected X,Y
