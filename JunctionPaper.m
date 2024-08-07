%% ADCP data: PAPER 3
clearvars
% Input parameters
lf = [500, 500, 50, 50, 500];   %Reg pars flow
ls = [100, 100];                %Reg pars salt
res = [50, 20];                 %Mesh resolution
evres = [24, 50, 20];          %Time x Lateral x Vertical - decomp resolution.

% Selection of data
dataset_idx = 2;
[adcp, ctd, h, dataset_name, ~] = import_data2(dataset_idx);
for ti = 1:3
    % Solve for flow and salinity using the ADCPTools
    [flow{ti}, sal{ti}] = solve_flow_salinity(adcp, ctd, ti, dataset_name, h, lf, ls, res);
    % Evaluate on regular t-y-sig grid and get Decomposition object.
    [U{ti}, S{ti}, D{ti}] = get_us(flow{ti}, sal{ti}, evres);
end
  figure
for ti = 1:3
  subplot(3,2,2*ti-1)
    hold on
    for t = 1:evres(1)
        plot(D{ti}.X.y, squeeze(U{ti}{1}(t, :, 3)), 'b')
        plot(D{ti}.X.y, squeeze(U{ti}{1}(t, :, 17)), 'r')
    end
    hold off
    subplot(3,2,2*ti)
        hold on
    for t = 1:evres(1)
        plot(D{ti}.X.y, squeeze(S{ti}{1}(t, :, 3)), 'b')
        plot(D{ti}.X.y, squeeze(S{ti}{1}(t, :, 17)), 'r')
    end
    hold off
end
M2T = flow.solver.model.periods(1,1)/(3600*24);
visres = [12, evres(2), 3];
t = linspace(0, M2T, visres(1)); % in days
n = linspace(min(flow.solver.mesh.n_left)+.5, max(flow.solver.mesh.n_right)-.5, visres(2));
sig = linspace(.2, .8, visres(3));




DF=D.decompose_function(U);
D.plot_components(DF, "velmap");
DFs=D.decompose_function(S);
D.plot_components(DFs, "salmap");

[SF, ~] = D.decompose_product(U,U);
figure;
bar(SF)
[SF, ~] = D.decompose_product(U,S);
figure;
bar(SF)

% Decomposition, visualization and rendering paper figures.
%D.plot_components(DF, "velmap")
%D.plot_components(DFs, "salmap")

%flow.evaluate(1, )

% [DF, AF, DFf, AFf] = D.decompose_function(u);


%flow.plot_solution({flow_model.all_names{[1]}}, representation = "Aphi", sol_idx = 1) % u,v,w
%sal.plot_solution({sal_model.all_names{[1]}}, representation = "Aphi", sol_idx = 1) % u,v,w
%clim([0,15])
%colormap(helpers.cmaps("salmap"))

%Dsums = sum(cat(4, DFfs{:}), 4); %time x lateral x vertical

%animate_solution(Dsum, X, nqt, Dsums);