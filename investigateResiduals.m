function investigateResiduals()

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