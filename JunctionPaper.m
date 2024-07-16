%% ADCP data: PAPER 3
clearvars
% Input parameters
lf = [500, 500, 50, 50, 500];   %Reg pars flow
ls = [100, 100];                %Reg pars flow
res = [50, 20];                 %Mesh resolution
evres = [20, 20, 20];          %Time x Lateral x Vertical - decomp resolution.

% Selection of data
dataset_idx = 1;
tak_idx = 1;
[adcp, ctd, h, dataset_name, ~] = import_data2(dataset_idx);
% Solve for flow and salinity using the ADCPTools
[flow, sal] = solve_flow_salinity(adcp, ctd, tak_idx, dataset_name, h, lf, ls, res);
[u, v, w, s, D] = get_us(flow, sal, evres);

% Decomposition, visualization and rendering paper figures.
D.plot_components(DF, "velmap")
D.plot_components(DFs, "salmap")
% [DF, AF, DFf, AFf] = D.decompose_function(u);


%flow.plot_solution({flow_model.all_names{[1]}}, representation = "Aphi", sol_idx = 1) % u,v,w
%sal.plot_solution({sal_model.all_names{[1]}}, representation = "Aphi", sol_idx = 1) % u,v,w
%clim([0,15])
%colormap(helpers.cmaps("salmap"))

%Dsums = sum(cat(4, DFfs{:}), 4); %time x lateral x vertical

%animate_solution(Dsum, X, nqt, Dsums);