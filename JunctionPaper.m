%% ADCP data: PAPER 3
clearvars

% Aggregate regularization parameters
l = [50];

% Datasets to analyze
didx = 1;
bidx = 1:3;
[lf, ls, ures, sres, evres, dname, bname, flow, salt, F, D, DF, AF, NF, US, adcp, ctd, h] = preallo(l);

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
    Tlim{di} = get_Tlim(flow{di}, bidx); % obtain
    for bi = bidx
        [F{di}{bi}{1}, F{di}{bi}{2}, D{di}{bi}] = get_us(flow{di}{bi}, salt{di}{bi}, evres, 1:numel(l), Tlim{di});
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

%% Full decomposition and correlation analysis
% figure
ri=1;
for di = didx
    for bi = bidx
            UU{di}{bi}{ri,1} = D{di}{bi}.get_prod_components_all(DF{di}{bi}{1}{ri, 1}, DF{di}{bi}{1}{ri, 1}); % Flow * Flow
            SS{di}{bi}{ri,1} = D{di}{bi}.get_prod_components_all(DF{di}{bi}{2}{ri, 1}, DF{di}{bi}{2}{ri, 1}); % Salt * Salt
            SF{di}{bi}{ri,1} = D{di}{bi}.get_prod_components_all(DF{di}{bi}{1}{ri, 1}, DF{di}{bi}{2}{ri, 1}); % Salt Flux
            %D{di}{bi}.plotUSF(UU{di}{bi}{ri,1}, 0:7)
            %D{di}{bi}.plotUSF(SS{di}{bi}{ri,1}, 0:7)
            %D{di}{bi}.plotUSF(SF{di}{bi}{ri,1}, 0:7)
    end
end
%%

%% 1D tide: Resolved vs unresolved.
fig=makefigure(8,8);
[resolved, unresolved, perc] = plot_1DT(fig, didx, bidx, D, bi, SF, ri, bname);
%fontsize(fig, 10, "points")
%% Idealized 2DV
fig=makefigure(8,8);
[resolved, unresolved, perc] = plot_2DV(fig, didx, bidx, SF, ri, D, bname, resolved, unresolved, perc);
%fontsize(fig, 10, "points")
%% Plot 1D transports
plot_1D_transports(didx, bidx, SF, ri, D, bname, UU);
%% Plot 1D transports Salt Flux
%fig=makefigure(16,8)
plot_1D_tidalTransport(didx, bidx, SF, ri, D, bname, UU);
%fontsize(fig, 10, "points")
%% Plot Net Salt Flux
%fig=makefigure(18,8)

plot_NetSaltFlux(didx, bidx, SF, bname)
%fontsize(fig, 10, "points")

%% Plot Net Balance
 %fig=makefigure(8,14)

plot_NetSaltBalance(didx, bidx, SF, bname)
%fontsize(fig, 10, "points")

%% Plot some components t-sigma (z) and y-t
% quick_inspect_plot(2, [1, 3], lf, D, DF, AF, c);
%% Plot sigma averaged HC
treal = HAK_ty_plot(D, AF);
%% Plot y-averaged NM (t - z plot)
%NM_tsig_plot(treal, D, AF);
%% Section 1b of the results: Instantaneous 3D dynamics
evaluate_strat(didx, bidx, lf, F, D, AF);
%% Section 2 of the results: Norms of flow and salinity, INSTANTANEOUS
NF = norms_bar_instantaneous_plot(D, didx, bidx, lf, NF);
%% Section 2 of the results: Norms of flow and salinity, RESIDUAL
%norms_bar_plot(didx, lf, bidx, NF);
%% 1D plots
csa_strat_plot(didx, bidx, F, D, AF);
%% Plot subtidal dynamics with top view
residual_top_plot(didx, bidx, flow, salt, F, D, AF);

%% Make only legends and colorbars for salinity

make_legends_colorbars()



%% save figs
% save_figs("./bin")

