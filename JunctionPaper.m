%% ADCP data: PAPER 3
clearvars

% Aggregate regularization parameters
l = [50];

% Datasets to analyze
didx = 1;
bidx = 1:3;
[lf, ls, ures, sres, evres, dname, bname, flow, salt, F, D, DF, AF, NF, SF, adcp, ctd, h] = preallo(l);

for di = didx
    [adcp{di}, ctd{di}, h{di}] = import_data2(dname{di});
    for bi = bidx
        [flow{di}{bi}, salt{di}{bi}] = solve_flow_salinity(dname{di}, bi, adcp{di}, ctd{di}, h{di}, lf, ls, ures, sres); % Solve for flow and salinity using the ADCPTools
    end
end

%% Evaluate solution and compute instantaneous flux.

% Order of variables:
% dataset_idx - branch_idx - var_idx (u, s, f) - reg_idx
% E.G. F{1}{1}{1}{2, 3}  : dataset 1, branch 1, variable 1,
% reg_pars 2, component 3 (which is w)
for di = didx
    for bi = bidx
        [F{di}{bi}{1}, F{di}{bi}{2}, D{di}{bi}] = get_us(flow{di}{bi}, salt{di}{bi}, evres, 1:numel(l));
        for ri = 1:size(lf, 1)
            for c = 1:3 %three components of the velocity to obtain a 3D flux
                F{di}{bi}{3}{ri, c} = F{di}{bi}{1}{ri, c}.*F{di}{bi}{2}{ri, 1};
            end
        end
    end
end

%% Decompose solution - only u, s and f
% NH{di}{bi} = D{di}{bi}.decompose_product_diag(DH{di}{bi}, DH{di}{bi});

for di = didx
    for bi = bidx
        for vi = 1:3 %Flow, Salt, Flux
            for ri = 1:size(lf, 1)
                for c = 1:size(F{di}{bi}{vi}, 2)
                    %                     disp(c)
                    [DF{di}{bi}{vi}{ri, c}, AF{di}{bi}{vi}{ri, c}] = D{di}{bi}.decompose_function(F{di}{bi}{vi}{ri, c}); % U-Flow, Salinity, Flux
                    NF{di}{bi}{vi}{ri, c} = D{di}{bi}.decompose_product_diag(DF{di}{bi}{vi}{ri, c}, DF{di}{bi}{vi}{ri, c}); % Norm of U, S, F
                end
            end
        end
    end
end

%% CSA fig
csa_figs(didx, bidx, D, flow, AF, bname) % (didx, bidx, D, flow, AF, tak_name)

%% Net salt flux components.
for di = didx
    for bi = bidx
        for ri = 1:size(lf, 1)
            SF{di}{bi}{ri,1} = D{di}{bi}.decompose_product_diag(DF{di}{bi}{1}{ri, 1}, DF{di}{bi}{2}{ri, 1}); % Salt Flux
        end
    end
end

%% Section 1b of the results: Instantaneous 3D dynamics
vlimbot = [.125, .375];
vlimsur = [.625, .875];

% close all
% Instantaneous Function Bottom / Surface
for di = didx
    for bi = bidx
        for vi = 1:3 %Flow, Salt, Flux
            for ri = 1:size(lf, 1)
                for c = 1:size(F{di}{bi}{vi}, 2)
                    IFB{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(F{di}{bi}{vi}{ri, c}, vlimbot); % U-Flow, Salinity, Flux instantaneous
                    IFS{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(F{di}{bi}{vi}{ri, c}, vlimsur); % U-Flow, Salinity, Flux instantaneous
                    
                    AIFB{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(AF{di}{bi}{vi}{ri, c}{5}, vlimbot);
                    AIFS{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(AF{di}{bi}{vi}{ri, c}{5}, vlimsur);
                end
            end
        end
    end
end

%% Section 2 of the results: Norms of flow and salinity, and their
% product
[n, a] = D{didx(1)}{bidx(1)}.get_names("u");
tit = ["Flow", "Salt", "Flux"];


for di = didx
    for bi = bidx
        [fig, TL] = prep_fig([20, 20], [numel(l), 4]);
        for ri = 1:size(lf, 1)
            for vi = 1:3 %flow, salt, flux
                for c = 1
                    [NF{di}{bi}{vi}{ri, c}(:,2), NF{di}{bi}{vi}{ri, c}(:,3), NF{di}{bi}{vi}{ri, c}(:,4), NF{di}{bi}{vi}{ri, c}(:,5)] = ...
                        scale_sort(NF{di}{bi}{vi}{ri, c}(:,1)); % normalized_sorted
                    % , idx, normalized_sorted_cumsum, normalized
                    nexttile;
                    bar(NF{di}{bi}{vi}{ri, c}(:,2))
                    xticklabels(a(NF{di}{bi}{vi}{ri, c}(:,3)))
                    title([tit(vi)])
                end

            end
            nexttile;
            bar(SF{di}{bi}{ri})
            %             xticklabels(a(NF{1,bi, di, ri}(:,3)))
            title("Salt Flux")
        end

        sgtitle(bname{di}{bi})
    end
end

%% Plot all components

cm = ["velmap", "salmap", "fluxmap"];

for di = didx
    for bi = bidx
        for ri = 1:size(lf, 1)
            for vi = 1:3 %flow, salt, flux
                D{di}{bi}.plot_components(AF{di}{bi}{vi}{ri, c}, cm(vi))
                sgtitle(["Averaged ", tit(vi), bname{di}{bi}, "with l = ", num2str(l(ri))])

                D{di}{bi}.plot_components(DF{di}{bi}{vi}{ri, c}, cm(vi))
                sgtitle(["Orthogonal ", tit(vi), bname{di}{bi}, "with l = ", num2str(l(ri))])
            end
        end
    end
end

%% Plot subtidal dynamics with top view

% Preliminary plot
rp = 1;
% One figure per dataset, one subplot per time
close all
yc = 1:9:evres(2); % coarsening of lateral profiles.

maxs{1} = 3; % max depth-averaged salinities -> integers
maxs{2} = 14;
figure;
for di = didx
    %     prep_fig([10,10], [3,2])
%     
        ax=nexttile;
        ti=1;
        for bi = bidx

            xb_ = [flow{di}{bi}.solver.mesh.x_middle(1), flow{di}{bi}.solver.mesh.x_middle(end)];
            yb_ = [flow{di}{bi}.solver.mesh.y_middle(1), flow{di}{bi}.solver.mesh.y_middle(end)];

            xq = linspace(xb_(1), xb_(2), D{di}{bi}.sz(2));
            yq = linspace(yb_(1), yb_(2), D{di}{bi}.sz(2));

            hold on
            var = 2; %salt
            s0 = AF{di}{bi}{var}{rp,1}{3}(ti, :);
            sb = AIFB{di}{bi}{var}{rp,1}(ti, :);
            ss = AIFS{di}{bi}{var}{rp,1}(ti, :);


            ssc = 20;
            ax = append_salt(ax, salt{di}{bi}, sb, xq, yq, ssc, maxs{di});
            ax = append_salt(ax, salt{di}{bi}, s0, xq, yq, ssc, maxs{di});
            ax = append_salt(ax, salt{di}{bi}, ss, xq, yq, ssc, maxs{di});
%             colorbar;


            var = 1; % flow
            % subtidal: index 3 in AF
            [u0, v0] = flow{di}{bi}.solver.xs.sn2xy_vel(AF{di}{bi}{var}{rp,1}{3}(ti, :), AF{di}{bi}{var}{rp, 2}{3}(ti, :)); % depth-averaged subtidal
            [ub, vb] = flow{di}{bi}.solver.xs.sn2xy_vel(AIFB{di}{bi}{var}{rp,1}(ti, :), AIFB{di}{bi}{var}{rp,2}(ti, :));
            [us, vs] = flow{di}{bi}.solver.xs.sn2xy_vel(AIFS{di}{bi}{var}{rp,1}(ti, :), AIFS{di}{bi}{var}{rp,2}(ti, :));
            % flow
            fsc = 30*ssc;
            quiver(ax, xq(yc), yq(yc), fsc*u0(yc), fsc*v0(yc), 'k', 'LineWidth', 2, 'AutoScale', 'off') % Subtidal flow
            quiver(ax, xq(yc+2), yq(yc+2), fsc*ub(yc+2), fsc*vb(yc+2), 'r', 'LineWidth', 2, 'AutoScale', 'off') % Subtidal flow
            quiver(ax, xq(yc+4), yq(yc+4), fsc*us(yc+4), fsc*vs(yc+4), 'b', 'LineWidth', 2, 'AutoScale', 'off') % Subtidal flow

            hold on
            plot(ax, xb_, yb_, "*k")
            plot(ax, xq, yq)
            hold off
        end
        axis square
    
end


%% Plot tidal dynamics with top view


% Preliminary plot
rp = 1;
% One figure per dataset, one subplot per time
close all
yc = 1:9:evres(2); % coarsening of lateral profiles.

tplot = 1:4:24; % indices, so integers >0 -> approx 2 hours between each plot.
maxs{1} = 6; % max depth-averaged salinities -> integers
maxs{2} = 18;

for di = didx
    %     prep_fig([10,10], [3,2])
    figure;
    for ti = tplot
        ax=nexttile;
        for bi = bidx

            xb_ = [flow{di}{bi}.solver.mesh.x_middle(1), flow{di}{bi}.solver.mesh.x_middle(end)];
            yb_ = [flow{di}{bi}.solver.mesh.y_middle(1), flow{di}{bi}.solver.mesh.y_middle(end)];

            xq = linspace(xb_(1), xb_(2), D{di}{bi}.sz(2));
            yq = linspace(yb_(1), yb_(2), D{di}{bi}.sz(2));

            hold on
            var = 2; %salt
            s0 = AF{di}{bi}{var}{rp,1}{7}(ti, :);
            sb = IFB{di}{bi}{var}{rp,1}(ti, :);
            ss = IFS{di}{bi}{var}{rp,1}(ti, :);


            ssc = 20;
            ax = append_salt(ax, salt{di}{bi}, sb, xq, yq, ssc, maxs{di});
            ax = append_salt(ax, salt{di}{bi}, s0, xq, yq, ssc, maxs{di});
            ax = append_salt(ax, salt{di}{bi}, ss, xq, yq, ssc, maxs{di});



            var = 1; % flow
            [u0, v0] = flow{di}{bi}.solver.xs.sn2xy_vel(AF{di}{bi}{var}{rp,1}{7}(ti, :), AF{di}{bi}{var}{rp, 2}{7}(ti, :)); % depth-averaged tidal
            [ub, vb] = flow{di}{bi}.solver.xs.sn2xy_vel(IFB{di}{bi}{var}{rp,1}(ti, :), IFB{di}{bi}{var}{rp, 2}(ti, :));
            [us, vs] = flow{di}{bi}.solver.xs.sn2xy_vel(IFS{di}{bi}{var}{rp,1}(ti, :), IFS{di}{bi}{var}{rp, 2}(ti, :));
            % flow
            fsc = 30*ssc;
            quiver(ax, xq(yc), yq(yc), fsc*u0(yc), fsc*v0(yc), 'k', 'LineWidth', 2, 'AutoScale', 'off') % Subtidal flow
            quiver(ax, xq(yc+2), yq(yc+2), fsc*ub(yc+2), fsc*vb(yc+2), 'r', 'LineWidth', 2, 'AutoScale', 'off') % Subtidal flow
            quiver(ax, xq(yc+4), yq(yc+4), fsc*us(yc+4), fsc*vs(yc+4), 'b', 'LineWidth', 2, 'AutoScale', 'off') % Subtidal flow


            hold on
            plot(ax, xb_, yb_, "*k")
%             plot(ax, xb2, yb2, "*k")
            plot(ax, xq, yq)
            hold off
        end
        axis square

    end

    %     legend(["u-bot", "u-surf", "u-avg", "s-bot", "s-sur", "s-avg"])
end