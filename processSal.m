function S = processSal(U, ctd, DS)

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

% plot(U.adcpI)
% hold on
% plot(S.i)
% disp(numel(unique(adcpH)))
% disp(numel(unique(S.t)))
% disp(max(adcpH) - max(S.t))
% disp(min(adcpH) - min(S.t))
% datestr(max(adcpH)- min(adcpH))

% t2 = datenum(adcpT);
% 
% figure;
% plot(adcpH)
% hold on
% plot(S.t)
S.eta = U.eta(S.t); % Waterlevel per CTD sample
% plot(S.eta)
% 
% figure;
% subplot(211)
% plot(ctd.X, ctd.Y)
% 
% subplot(212)
% try
%     plot(templon, templat)
% 
% figure
% plot(S.X, S.Y, '*')
% hold on
% plot(xnew(:,1), xnew(:,2), '*')
% U.B.plot()
% U.xs.plot()
% catch
%     
% end

S.Z          =  - ctd.Z + S.eta;



% 
% temp = datetime(S.t, 'ConvertFrom','datenum','Format','HH:mm:ss');
% 
% [S.hr, ind] = getHour(S.t, hour/2);

[S.s, S.n] = U.xs.xy2sn(S.X, S.Y);

S.zb = n2zb(S.n, U.mesh_mean);
S.sig = z2sig(S.Z, S.zb, S.eta); % IS THIS ALLOWED IN SIGMA COORDINATES?? - TAKE CARE
S.Z = sig2z(S.sig, S.zb, S.eta);

% figure
% plot(S.sig)
% figure;
% % Inspect CTD data on mesh
% plot(S.n, S.Z,'*')
% hold on
% U.mesh_mean.plot()



[S.Xp, S.Yp] = U.xs.sn2xy(0*S.s, S.n); % Project on cross section by setting s = 0. [Xp, Yp] = projected X,Y


% Now do the actual, hour-based calculation.



% 
% [S.sp, S.np] = U.xs.xy2sn(S.Xp, S.Yp);
% 
% S.sn = [S.t, S.s, S.n, S.Z];
% S.snp = [S.t, S.sp, S.np, S.Z];
% 
% S.snN = rescale(S.sn, 'InputMin', min(S.sn), 'InputMax', max(S.sn));
% S.snpN = rescale(S.snp, 'InputMin', min(S.snp), 'InputMax', max(S.snp));
% 
% figure;
% scatter3(S.snN(:,1), S.snN(:,3), S.snN(:,4))
% axis equal
% xlabel('t'); ylabel('n'); zlabel('z')
% 
% figure;
% scatter3(S.snpN(:,1), S.snpN(:,3), S.snpN(:,4))
% axis equal
% xlabel('t'); ylabel('n'); zlabel('z')

% SI = LoessInterpolator();
% SI.known = [S.snpN S.S]';
% SI.span = .1;
% SI.order = 1;
% SI.robust_iterations = 1;
% 
% nl = [0, logspace(-2, 2, 29)];
% nnl = length(nl);
% nums = numel(S.S);
% 
% for i = 1:10
%     epsi = randn(size(S.snpN));
%     disp(nl(i))
%     tic
%     SI1(:,i) = SI.interpolate(S.snpN' + nl(i)*epsi');
%     t1 = toc;
%     tic
%     SI2(:,i) = griddatan(S.snpN, S.S, S.snpN + nl(i)*epsi); %150 times slower than loess in standard config
%     t2 = toc;
%     nor(i,1) = i;
% 
%     nor(i,2) = sum(abs(S.S - SI1(:,i)).^2, 'omitnan');
%     nor(i,3) = sum(isnan(SI1(:,i)))/nums;
%     nor(i,4) = t1;
% 
%     nor(i,5) = sum(abs(S.S - SI2(:,i)).^2, 'omitnan');
%     nor(i,6) = sum(isnan(SI2(:,i)))/nums;
%     nor(i,7) = t2;
%     disp(['LOESS and GDN comparison: ', num2str(nor(i,:))])
% end
% 
% 
% [SS, ind] = sort(S.S);
% plot(SS)
% hold on
% plot(SI1(ind,i))
% plot(SI2(ind,i))
% legend('Salt', 'SI1', 'SI2')


% [x, y] = meshgrid(1:nums, 1:nnl);
% subplot(211)
% surf(x,y, abs(repmat(S.S, 1, nnl) - SI1), '--')
% subplot(212)
% surf(x,y, abs(repmat(S.S, 1, nnl) - SI2), '--')
% % plot(S.S)

% figure;
% 
% scatter3(S.t, S.s, S.n)
% xlabel('t'); ylabel('s'); zlabel('n');

% figure;
% subplot(221)
% scatter3(S.t, S.s, S.n)
% xlabel('t'); ylabel('s'); zlabel('n');
% 
% subplot(222)
% scatter3(S.t, 0*S.s, S.n)
% xlabel('t'); ylabel('0'); zlabel('n');
% 
% subplot(223)
% scatter3(S.t, S.X, S.Y)
% xlabel('t'); ylabel('X'); zlabel('Y');
% 
% subplot(224)
% scatter3(S.t, S.Xp, S.Yp)
% xlabel('t'); ylabel('Xp'); zlabel('Yp');
% S.P_tmp{hour}      = [S.x{hour},S.y{hour},S.Z{hour}];

% create CTD position data vector
%     cur_mesh = U.mesh(hour);
%
%     % Transform the CTD data (UTM) to mesh coordinates. Firstly, transform
%     % to the same coordinates as the mesh cells - what are they?
%     %P          =
%     % [ss,ns] = xs.xy2sn(cur_mesh.x_middle(cur_mesh.col_to_cell),cur_mesh.y_middle(cur_mesh.col_to_cell),orig_vel{crp}(:,1),orig_vel{crp}(:,2));
%     [CTD.sctd{hour}, CTD.nctd{hour}] = U.xs.xy2sn( CTD.x{hour}, CTD.y{hour});
%     CTD.zctd{hour} = CTD.Z{hour};
%
%     CTD.mn{hour} = mean(CTD.nctd{hour});
%     CTD.ms{hour} = mean(CTD.sctd{hour});
%     CTD.velh{hour} = U.vel{hour};
%     CTD.n_cells{hour} = mean(cur_mesh.n_patch,1)';
%     CTD.z_cells{hour} = mean(cur_mesh.z_patch,1)';
%
%     CTD.Fv{hour}  = scatteredInterpolant(CTD.sctd{hour} - CTD.ms{hour}, CTD.nctd{hour} - CTD.mn{hour}, CTD.zctd{hour}, CTD.sal{hour}, 'linear', 'nearest');
%     CTD.S{hour} = CTD.Fv{hour}(0.*ones(size(CTD.n_cells{hour})),  CTD.n_cells{hour},  CTD.z_cells{hour}); %here, take other coordinates

%eta(hour) = mesh(hour).water_level;
end