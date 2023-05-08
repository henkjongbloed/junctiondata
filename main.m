%% ADCP data: Method paper
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
close all % Close all figures
figsave = 0;

RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

ds = 2;
tak = 1;

DS = {'OMHA14','NMOMNW15'}; DS = DS{ds};
BN = {{' Hartel Canal'; ' Old Meuse South'; ' Old Meuse North'},...
    {' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}}; BN = BN{ds};

BNs = {{'HK2'; 'OMS2'; 'OMN2'},...
    {'NM2'; 'OM2'; 'RW2'}}; BNs = BNs{ds};

%% Import preprocessed semi-raw data (ADCP / CTD / H structs)

[adcp, ctd, h] = import_data(RF, DS);

V = rdi.VMADCP(adcp{tak}); % create VMADCP object
% V.horizontal_position_provider(1) = [];
mtc= MagneticDeviationTwoCycle;
mtc.plot_tracks(V)
mtc.estimate_deviation(V)
mtc.estimate_deviation(V,true)

constituents = {'M2', 'M4'};
% 
water_level = VaryingWaterLevel(datetime(datevec(h.t)), h.wl); % set water level to varying
water_level.model = TidalModel(constituents = constituents);
water_level.model.components = {'eta'}; % Scalar
water_level.get_parameters();

V.water_level_object = water_level;
%% create bathymetry
bathy = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)

dat_selection = "default";
switch dat_selection
    case "manual"
        [ef, xs] = cross_section_selector(V);
    case "default"
        xs = XSection(V); % define the cross-section
        ef = EnsembleFilter(V);
    case "saved" % Save ensemblefilter and cross-section once
        xs = XSection(V); % define the cross-section
        ef = EnsembleFilter(V);
    otherwise 
        xs = XSection(V); % define the cross-section
        ef = EnsembleFilter(V);
end

if any(strcmp({' New Meuse'; ' Rotterdam Waterway';' Hartel Canal'}, BN))
    xs.revert();
end

U.BN = BN;

if any(strcmp({' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}, BN))
    U.DS = 2015;
else
    U.DS = 2014;
end


%% create mesh

mesh_maker = SigmaZetaMeshFromVMADCP(V, bathy, xs); % mesh generator
mesh_maker.deltan = 100;
mesh_maker.deltaz = 2;
mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level
Bw = mesh_mean.nw(2) - mesh_mean.nw(1); % Little fishy but it works
Hw = abs(min(mesh_mean.zb_all));
mesh_maker.deltan = (Bw + 1)/21;
mesh_maker.deltaz = (Hw+1)/11;
mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level
% mesh_mean.plot3()
% hold on
% bathy.plot


%tid = TidalModel(constituents = {'M2', 'M4'});
%tay = TaylorModel(s_order = [1,1,1], n_order = [1,1,1], sigma_order = [1,1,1]);

%data_model = TaylorTidalModel(taylor = tay, tidal = tid, rotation = pi/2);

% data_model = TaylorTidalModel(constituents = constituents, s_order = [1,1,1], n_order = [1,1,1], sigma_order = [1,1,1]);
%phi = linspace(-pi, pi, 5);
% Benchmark: classic fit
% classic_model = TidalModel(constituents = constituents, rotation = 0);
% 
% % classic_model.rotation = -phi;
% T = LocationBasedVelocitySolver(V, mesh_mean, bathy, xs, ef, classic_model);
% pars0 = T.get_parameters();
% figure;
% tiledlayout('flow');nexttile;
% for i = 1:size(pars0, 2)
%     mesh_mean.plot(pars0(:,i));colorbar;nexttile;
% end
% 
% 
% % Benchmark: classic fit
% classic_model = TidalModel(constituents = constituents, rotation = phi);
% 
% % classic_model.rotation = -phi;
% T = LocationBasedVelocitySolver(V, mesh_mean, bathy, xs, ef, classic_model);
% pars1 = T.get_parameters();
% figure;
% tiledlayout('flow');nexttile;
% for i = 1:size(pars1, 2)
%     mesh_mean.plot(pars1(:,i));colorbar;nexttile;
% end
% 
% 
% % Benchmark: classic fit
% classic_model = TidalModel(constituents = constituents, rotation = -phi);
% 
% % classic_model.rotation = -phi;
% T = LocationBasedVelocitySolver(V, mesh_mean, bathy, xs, ef, classic_model);
% pars2 = T.get_parameters();
% figure;
% tiledlayout('flow');nexttile;
% for i = 1:size(pars2, 2)
%     mesh_mean.plot(pars2(:,i));colorbar;nexttile;
% end
data_model = TaylorTidalModel(constituents = constituents, s_order = [1,1,1], n_order = [1,1,1], sigma_order = [1,1,1]);

opts = SolverOptions();
opts.algorithm = "pcg";
opts.pcg_iter = 4000;
opts.pcg_tol = 1e-10;
opts.reg_pars = {0*[1,1,1,1,1], 1*[100,100,10,10,100]};

R = Regularization(bathy = bathy, xs = xs, mesh = mesh_mean, model = data_model);
R.assemble_matrices()

data_model.rotation = xs.angle();
S = LocationBasedVelocitySolver(V, mesh_mean, bathy, xs, ef, data_model, R, opts);
MP = S.get_parameters();

%mesh_mean.plot(MP.pars(:,3))
%colorbar;


% MP.cross_validate();

figure;
MP.plot_solution({data_model.names{1}{1:20}}, 1:2);

figure;
MP.plot_solution({data_model.names{2}{1:20}}, 1:2);

figure;
MP.plot_solution({data_model.names{3}{1:20}}, 1:2);


%%

%[adcp, ctd, h] = import_data_(RF, DS);
%adcp{1} = raw_dat;
% %% Data analysis
% plotvars = 1;
% for j = 1:length(tak)
%     b = tak(j);
%     [V{j}, U{j}] = processVel(adcp{b}, h, BN{b}); % V = VMADCP objects, U = generic velocity data (plus more)
%     
%     %     S{j} = processSal(U{j}, ctd{b}, DS); %S = Salinity data
%     %     F{j} = processHours(U{j}, S{j});
%     %     F{j} = processFlux(F{j});
%     %     [U{j}.I_loc, U{j}.I] = get_MoransI(U{j}, dat.IM, 1);
%     %     [U{j}.I_loc0, U{j}.I0] = get_MoransI(U{j}, dat.IM, 0);
%     %     [U{j}.I_locD, U{j}.ID] = get_MoransI(U{j}, dat.IM, -1);
% 
%     %
%     reg_dat{j} = dat;
% 
%     est_opts.generate = 'local';
%     est_opts.noise_levels = linspace(0, 1, 12);
% %     reg_dat{j} = est_parameters(reg_dat{j}, est_opts);
%     %     reg_dat{j}.names =
%     plt = 1;
% 
%     if plt
%         noise_ind = 1:12;
% %                  plotPe(reg_dat{j})
% %         makefigure(21,29.7/2)
% % 
% %         plotSe(reg_dat{j}, noise_ind, U{j}.T.velocity_model.names, est_opts)
% %         makefigure(21,29.7/2)
% %         ld_idx = 1:12;%[1, 4, 8, 12];
% %         plotNoiseLC(reg_dat{j}, ld_idx, est_opts)
% %         makefigure(21,29.7/2)
% %         lc_idx = 1:12;
% %         plotNoiseLD(reg_dat{j}, lc_idx, est_opts)
%           tplot = 0:360:13*3600;
%           %u = p2u(dat.p0, U{1}.T.velocity_model.names, U{1,1}.T.adcp.water_level_object, U{1,1}.mesh_mean, tplot);
% %                  plotSe(reg_dat{j}, noise_ind, U{j}.T.velocity_model.names, est_opts)
% 
% 
%         %         plotddp(dat.p0,U{1}.T.velocity_model.names([1:5, 16:20]), U{1}.T.velocity_model.names,  U{1}.mesh_mean, 1, U{1}.velmap)
%         %         plotp(dat.p(:,1,1), U{1}.T.velocity_model.names([1:5, 16:20]), U{1}.T.velocity_model.names,  U{1}.mesh_mean, 1, U{1}.velmap)
%         %
%         %
% %         plotp(dat.p0 ,U{1}.T.velocity_model.names, U{1}.T.velocity_model.names, 6, U{1})
%         %         plotp(dat.p(:,end,end),U{1}.T.velocity_model.names, U{1}.T.velocity_model.names,  U{1}.mesh_mean, 6, U{1}.velmap)
%         %        clear p_cell;
% %         plotp(p_cell{1,1}, U{1}.T.velocity_model.names, U{1}.T.velocity_model.names,  U{1}.mesh_mean, 6, U{1}.velmap)
%         %         plotp(pab2pAphi(dat.p(:,end,end),U{1}.T.velocity_model.names), U{1}.T.velocity_model.names, U{1}.T.velocity_model.names,  U{1}.mesh_mean, 6, U{1}.velmap)
% 
%         %           [p_cell{1,2},~] = pab2pAphi(dat.p(:,19,1), U{1}.T.velocity_model.names);
% %         [p_cell{1,1}, tid_names] = pab2pAphi(mean(squeeze(dat.p(:,1,:)),2), U{1}.T.velocity_model.names);
% %         [p_cell{1,2}, ~] = pab2pAphi(mean(squeeze(dat.p(:,6,:)),2), U{1}.T.velocity_model.names);
% %         [p_cell{1,3}, ~] = pab2pAphi(mean(squeeze(dat.p(:,9,:)),2), U{1}.T.velocity_model.names);
% % %         makefigure(21,29.7/2)
%         names = U{1}.T.velocity_model.names;
% %             p0 = pab2pAphi(dat.p0, U{1}.T.velocity_model.names);
% %         makefigure(50, 22)
% %          plotp(dat.p1, names, names,  10, U{1})
%         [p1(:,1), tid_names] = pab2pAphi(dat.p1(:,1), U{1}.T.velocity_model.names);
%         p1(:,2) = pab2pAphi(dat.p1(:,2), U{1}.T.velocity_model.names);
%         p1(:,3) = pab2pAphi(dat.p1(:,3), U{1}.T.velocity_model.names);
% 
%         plotp_compare({p1(:,1), p1(:,2), p1(:,3)},tid_names([1:4, 5, 9]), tid_names,  U{1})
% 
%         %         makefigure(21,29.7/2)
% %         plotp(p_cell{1,2},...
% %             tid_names, tid_names,  6, U{1})
% 
%         %([1, 2, 4, 21, 22, 24, 41, 42, 44])
% 
% %         plotp(mean(squeeze(dat.p(:,18,1:end)),2),U{1}.T.velocity_model.names, U{1}.T.velocity_model.names,  U{1}.mesh_mean, 10, U{1}.velmap)
%     end
% 
%     %     plotp(U{j}.mesh_mean, U{j}.T.velocity_model.names, pgen, U{j}.velmap)
%     %     plotp(U{j}.mesh_mean, U{j}.T.velocity_model.names, pest, U{j}.velmap)
% 
% 
%     %     res = reg_dat{j}.b - reg_dat{j}.M*reg_dat{j}.p; % Residuals
%     %     if plotvars
%     %         figure;
%     %         subplot(2,2,[1,3])
%     %         hold on
%     %         grid on
%     %         view(30,40)
%     %
%     %         for c = 1:U{j}.mesh_mean.ncells % residuals per mesh cell
%     %             resCell{c,1} = res(c == reg_dat{j}.cell_idx);
%     %             resCell_mean(c,1) = mean(resCell{c,1});
%     %             resCell_std(c,1) = std(resCell{c,1});
%     %             nbins = 20;
%     %
%     %             [res_N(c,:), res_edges(c,:)] = histcounts(resCell{c,1},nbins, 'Normalization', 'pdf');
%     %             dr = res_edges(c,2) - res_edges(c,1);
%     %             dr_plot = res_edges(c,1:end-1) + dr;
%     %
%     %             %         plot3(dr_plot, c*ones(length(dr_plot)), res_N(c,:)) % plot in 3D space
%     %             dr_plot2 = linspace(min(dr_plot), max(dr_plot), 100);
%     %             plot3(dr_plot2, c*ones(length(dr_plot2)), normpdf(dr_plot2, resCell_mean(c,1), resCell_std(c,1)))
%     %             plot3(resCell_mean(c,1), c, 0, '*', 'MarkerSize', 3)
%     %             %         errorbar(x,y,err,'horizontal')
%     %         end
%     %         res_mean = mean(res);
%     %         res_std = std(res);
%     %         [res_N(c+1,:), res_edges(c+1,:)] = histcounts(res, nbins, 'Normalization', 'pdf');
%     %         dr = res_edges(c+1,2) - res_edges(c+1,1);
%     %         dr_plot = res_edges(c+1,1:end-1) + dr;
%     %
%     %         %         plot3(dr_plot, c*ones(length(dr_plot)), res_N(c,:)) % plot in 3D space
%     %         dr_plot2 = linspace(min(dr_plot), max(dr_plot), 100);
%     %         plot3(dr_plot2, 0*ones(length(dr_plot2)), normpdf(dr_plot2, res_mean, res_std), 'LineWidth', 3)
%     %         plot3(res_mean, 0, 0, '*', 'MarkerSize', 5)
%     %         % Plot std ranges
%     %
%     %         plot3([res_mean + res_std; resCell_mean + resCell_std], 0:c, 0*[0:c])
%     %         plot3([res_mean - res_std; resCell_mean - resCell_std], 0:c, 0*[0:c])
%     %         hold off
%     %
%     %         %     figure;
%     %         subplot(2,2,2)
%     %         U{j}.mesh_mean.plot(resCell_mean);
%     %         amax = max(abs(resCell_mean), [], 'omitnan');
%     %         caxis([-amax, amax])
%     %         colormap(U{j}.velmap)
%     %         set(gca, 'XDir','reverse') % Very important
%     %         title('mean residuals')
%     %         colorbar;
%     %         subplot(2,2,4)
%     %         colormap(U{j}.velmap(10:end, :))
%     %         U{j}.mesh_mean.plot(resCell_std);
%     %         %     amax = max(abs(resCell_std), [], 'omitnan');
%     %         %     caxis([0, amax])
%     %         set(gca, 'XDir','reverse') % Very important
%     %         title('std residuals')
%     %         colorbar;
%     %
%     %
%     %         figure;
%     %         hold on
%     %         plot(U{j}.I)
%     %         plot(U{j}.I0)
%     %         plot(U{j}.ID)
%     %
%     %         title([BN{j},': Morans-I score autocorrelation'])
%     %         xlabel('Parameter'); ylabel('Morans-I score');
%     %         legend(['Regularized: ', num2str(mean(U{j}.I))],...
%     %             ['Data only', num2str(mean(U{j}.I0))], ['Difference', num2str(mean(U{j}.ID))])
%     %     end
%     figsave = 0;
%     if figsave
%         figs =  findobj('type','figure');
%         nfig = length(figs);
%         FF = 'C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\figs\dec22\'; %RootFolder
%         %cd([RF,'WP2\TwoJunctions\code'])
%         % cd FF
%         for fi = 1:nfig
%             %             f = gcf;
%             %             st = f.Children.Title.String;
%             dest = fullfile(FF,[BNs{tak(j)}, num2str(fi),  '.png']);
%             %         set(gcf, 'Position', get(0, 'Screensize'));
%             saveas(figs(fi), dest);
%         end
%         close all
%     end
% end
% 
% 
% 
% %plot(U{1, 1}.T.velocity_model.names, U{1,1}.I)
% 
% 
% % [M1, b1, Mj1, bj1] = GD2M(gather_dat);
% 
% % condest(M1)
% % condest(M1'*M1)
% 
% % for j = 1:length(Mj)
% %     frm(j) = norm(Mj{j} - Mj1{j}, 'fro');
% %     frb(j) = norm(bj{j} - bj1{j}, 'fro');
% % end
% % plot_geom_all(U)
% %% Plotting
% % plt=1;
% % if plt
% %     tak=1;
% %     %     plotCovParsx(U{i})
% %     plotPars_all(U{j}, {'Old','New'})
% % end
% %%
% % figure;
% % plot(cve)
% % figure;
% % plot(p)
% 
% %
% % plotMap(U,S);
% %
% % plotDA_Flux(tak, U, F);
% % plotSubtidalFlow(tak, U);
% % plotSubtidalSalt(tak, F, U);
% % plot_geom_all(U, S)
% % % plotPars(tak, U);
% % plotT = 1:2:13;
% % plotTimeUS(tak, plotT, F, U);
% % plotCSA_Flux(tak, F, U);
% % plotWL_CSA(F, tak, BN);
% % %
% % % %% Tidal Ellipses
% % % %
% % plotTidalEllipses(U, F)
% % plotTidalArrows(U, F)
% %
% % % plotTZ(U,F, -.6)
% % plotTZ(U,F, 0)
% % % plotTZ(U,F, .6)
% %
% % Save all figs automatically
% 
% %
% %
