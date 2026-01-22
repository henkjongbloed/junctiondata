%% ADCP data: Hartel canal
clearvars

% Aggregate regularization parameters
l = [50];

% Datasets to analyze
didx = 1;
bidx = 1;
[lf, ls, ures, sres, evres, dname, bname, flow, salt, F, D, DF, AF, NF, US, adcp, ctd, h] = preallo(l);
 ures = [20,7];
for di = didx
    [adcp{di}, ctd{di}, h{di}] = import_data2(dname{di});
    for bi = bidx
        [flow{di}{bi}, salt{di}{bi}] = solve_flow_salinity(dname{di}, bi, adcp{di}, ctd{di}, h{di}, lf, ls, ures, sres); % Solve for flow and salinity using the ADCPTools
    end
end


sol = flow{1}{1};

%% From here, cross-validation (Ariane, start here)


N = [8,8]; % number of reg. pars to iterate over (N_cont, N_smooth)
K = 2; %k-fold cross-validation
sol.solver.opts.cv_iter=1;
[CV, rpc, rps] = sol.cross_validate_2D([0,0], [1e7, 1e7], N, K);

figure;
% make use of the ordering of the constraints: first continuity

contourf(helpers.symlog(rpc), helpers.symlog(rps), reshape([CV{:,1}]./CV{1,1}, N(2), N(1)), 100)
colorbar
colormap(flipud(helpers.cmaps('velmap')))
% clim([0,2]) % 0 - very good (too good to be true) 1
xlabel('cont lambda (symlog10)')
ylabel('smoothness lambda (symlog10)')
title('2D cross-validation: generalization error')