%% Read, filter, plot ADCP data
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
close all % Close all figures
figsave = 0;
%clc
RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

% ds = 2; %Choose which dataset to analyze
% ds = input('Enter dataset'); % 1 XOR 2
% tak = input('Enter branch'); % 1 or 1:2 or 2:3 or 1:3 (all)

ds = 1;
tak = 1;

DS = {'OMHA14','NMOMNW15'}; DS = DS{ds};
BN = {{' Hartel Canal'; ' Old Meuse South'; ' Old Meuse North'},...
    {' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}}; BN = BN{ds};

%% Import preprocessed semi-raw data (ADCP / CTD / H structs)

[adcp, ctd, figs] = import_data(RF, DS);

%% Data analysis
for j = 1:length(tak)
    b = tak(j);
    [V{j}, U{j}] = processVel(adcp{b}, figs, BN{b}); % V = VMADCP objects, U = generic velocity data (plus more)
    %     S{j} = processSal(U{j}, ctd{b}, DS); %S = Salinity data
    %     F{j} = processHours(U{j}, S{j});
    %     F{j} = processFlux(F{j});
    [U{j}.I_loc, U{j}.I] = get_MoransI(U{j}, dat.IM, 1);
    [U{j}.I_loc0, U{j}.I0] = get_MoransI(U{j}, dat.IM, 0);
    [U{j}.I_locD, U{j}.ID] = get_MoransI(U{j}, dat.IM, -1);


    reg_dat{j} = dat;
    
    res = reg_dat{j}.b - reg_dat{j}.M*reg_dat{j}.p; % Residuals
    figure;
    subplot(2,2,[1,3])
    hold on
    grid on
    view(30,40)

    for c = 1:U{j}.mesh_mean.ncells % residuals per mesh cell
        resCell{c,1} = res(c == reg_dat{j}.cell_idx);
        resCell_mean(c,1) = mean(resCell{c,1});
        resCell_std(c,1) = std(resCell{c,1});
        nbins = 20;

        [res_N(c,:), res_edges(c,:)] = histcounts(resCell{c,1},nbins, 'Normalization', 'pdf');
        dr = res_edges(c,2) - res_edges(c,1);
        dr_plot = res_edges(c,1:end-1) + dr;

%         plot3(dr_plot, c*ones(length(dr_plot)), res_N(c,:)) % plot in 3D space
        dr_plot2 = linspace(min(dr_plot), max(dr_plot), 100);
        plot3(dr_plot2, c*ones(length(dr_plot2)), normpdf(dr_plot2, resCell_mean(c,1), resCell_std(c,1)))
        plot3(resCell_mean(c,1), c, 0, '*', 'MarkerSize', 3)
%         errorbar(x,y,err,'horizontal')
    end
    res_mean = mean(res);
    res_std = std(res);
    [res_N(c+1,:), res_edges(c+1,:)] = histcounts(res, nbins, 'Normalization', 'pdf');
    dr = res_edges(c+1,2) - res_edges(c+1,1);
    dr_plot = res_edges(c+1,1:end-1) + dr;

    %         plot3(dr_plot, c*ones(length(dr_plot)), res_N(c,:)) % plot in 3D space
    dr_plot2 = linspace(min(dr_plot), max(dr_plot), 100);
    plot3(dr_plot2, 0*ones(length(dr_plot2)), normpdf(dr_plot2, res_mean, res_std), 'LineWidth', 3)
    plot3(res_mean, 0, 0, '*', 'MarkerSize', 5)
    % Plot std ranges

    plot3([res_mean + res_std; resCell_mean + resCell_std], 0:c, 0*[0:c])
    plot3([res_mean - res_std; resCell_mean - resCell_std], 0:c, 0*[0:c])
    hold off

%     figure;
    subplot(2,2,2)
    U{j}.mesh_mean.plot(resCell_mean);
        amax = max(abs(resCell_mean), [], 'omitnan');
        caxis([-amax, amax])
    colormap(U{j}.velmap)
    set(gca, 'XDir','reverse') % Very important
    title('mean residuals')
    colorbar;
    subplot(2,2,4)
    colormap(U{j}.velmap(10:end, :))
    U{j}.mesh_mean.plot(resCell_std);
%     amax = max(abs(resCell_std), [], 'omitnan');
%     caxis([0, amax])
    set(gca, 'XDir','reverse') % Very important
    title('std residuals')
    colorbar;


    figure;
    hold on
    plot(U{j}.I)
    plot(U{j}.I0)
    plot(U{j}.ID)

    title([BN{j},': Morans-I score autocorrelation'])
    xlabel('Parameter'); ylabel('Morans-I score');
    legend(['Regularized: ', num2str(mean(U{j}.I))],...
        ['Data only', num2str(mean(U{j}.I0))], ['Difference', num2str(mean(U{j}.ID))])
end



%plot(U{1, 1}.T.velocity_model.names, U{1,1}.I)


% [M1, b1, Mj1, bj1] = GD2M(gather_dat);

% condest(M1)
% condest(M1'*M1)

% for j = 1:length(Mj)
%     frm(j) = norm(Mj{j} - Mj1{j}, 'fro');
%     frb(j) = norm(bj{j} - bj1{j}, 'fro');
% end
% plot_geom_all(U)
%% Plotting
plt=1;
if plt
    tak=1;
    %     plotCovParsx(U{i})
    plotPars_all(U{j}, {'New'})
end
%%
% figure;
% plot(cve)
% figure;
% plot(p)

%
% plotMap(U,S);
%
% plotDA_Flux(tak, U, F);
% plotSubtidalFlow(tak, U);
% plotSubtidalSalt(tak, F, U);
% plot_geom_all(U, S)
% % plotPars(tak, U);
% plotT = 1:2:13;
% plotTimeUS(tak, plotT, F, U);
% plotCSA_Flux(tak, F, U);
% plotWL_CSA(F, tak, BN);
% %
% % %% Tidal Ellipses
% % %
% plotTidalEllipses(U, F)
% plotTidalArrows(U, F)
%
% % plotTZ(U,F, -.6)
% plotTZ(U,F, 0)
% % plotTZ(U,F, .6)
%
% Save all figs automatically
figsave = 0;
if figsave
    figs =  findobj('type','figure');
    nfig = length(figs);
    FF = 'C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\figs\21-10\'; %RootFolder
    %cd([RF,'WP2\TwoJunctions\code'])
    % cd FF
    for j = 1:nfig
        f = gcf;
        st = f.Children.Title.String;
        dest = fullfile(FF,[num2str(j), '.png']);
        %         set(gcf, 'Position', get(0, 'Screensize'));
        saveas(f, dest);
        close;
    end
end
%
%
