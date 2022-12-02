function dat = est_parameters(dat, est_opts)

M = dat.M; C = dat.C; Cg = dat.Cg; D = dat.D; D2 = dat.D2;
Np = size(M,2); % number of parameters
nc = max(dat.cell_idx); % number of cells
np = Np/nc; %number of parameters per cell
opts = dat.opts;
nepochs = opts.cv_iter; %nepochs now serves the role

if strcmp(est_opts.generate, 'nullspace')
    D0 = D;
    for row = 1:size(D0,1) % rescale smoothing matrix back - is this necessary?
        D0(row,:) = D0(row,:)/D(row,row);
    end
    tic;
    NS = generate_null_intersection({C, Cg, D0, D2});
    to = toc;
    fprintf('Finished calculating intersection of null spaces after %2.2f s \n', to)

    p = NS*randn(size(NS,2),nepochs);
elseif strcmp(est_opts.generate, 'local')
    p = repmat(dat.p0, 1, nepochs);
end
regP = combine_regpars(opts);
% p = nan([np,size(regP,1),nepochs]);

stdn = est_opts.noise_levels;
Mp = M'*M; Cp = C'*C; Cgp = Cg'*Cg; Dp = D'*D; D2p = D2'*D2;

B = M*p; % Unperturbed data
niter = length(stdn)*size(regP,1)*nepochs;
disp(niter)
i=0;
avg_err = nan([size(regP,1), length(stdn)]);
avg_perr = nan([np, size(regP,1), length(stdn)]);
pcg_opts = struct('michol','on','type','ict','droptol',1e-3);
for nn = 1:length(stdn)
    Bp = B + stdn(nn)*randn(size(B)); % Perturb measurements
    for rp = 1:size(regP,1)
        A = Mp + regP(rp,1)*Cp + regP(rp,2)*Cgp + regP(rp,3)*Dp + regP(rp,4)*D2p;
        pcg_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
        L = ichol(A, pcg_opts);
        for ep = 1:nepochs
            i = i+1;
            fprintf('Sensitivity analysis: %2.2f percent \n', 100*i/niter)
            phat(:, rp, ep) = pcg(A, M'*Bp(:,ep), 1e-8, size(A,2), L, L', p(:,ep)); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
            err(rp, ep) = sqrt(sum((p(:,ep)-phat(:, rp, ep) ).^2));
            rel_err(rp, ep) = err(rp, ep)./sqrt(sum(p(:,ep).^2));

            % Investigate sensitivity wrt all parameters small quantities
            for ip = 1:np
                perr(ip,rp, ep) = sqrt(sum((p(ip:np:end,ep)-phat(ip:np:end, rp, ep) ).^2));
                rel_perr(ip,rp, ep) = perr(ip,rp, ep)./sqrt(sum(p(ip:np:end,ep).^2));
            end
        end
        avg_err(rp, nn) = mean(err(rp, :));
        avg_perr(:, rp, nn) = mean(squeeze(perr(:, rp, :)),2);
        avg_rel_err(rp, nn) = mean(rel_err(rp, :));
        avg_rel_perr(:, rp, nn) = mean(squeeze(rel_perr(:, rp, :)),2);
    end
end

dat.avg_err = avg_err;
dat.avg_perr = avg_perr;
dat.avg_rel_err = avg_rel_err;
dat.avg_rel_perr = avg_rel_perr;
end