function RP = combine_regpars(opts)

% reg_pars must be a cell array of column vectors of not necessarily same size

if strcmp(opts.reg_vary, 'coupled')
    L = opts.reg_pars(1, [1,3]); % Only vary two reg. parameters
    n = length(L);
    [L{:}] = ndgrid(L{end:-1:1});
    L = cat(n+1,L{:});
    RP0 = fliplr(reshape(L,[],n));
    RP = [RP0 RP0];
    RP = RP(:, [1, 3, 2, 4]);
elseif strcmp(opts.reg_vary, 'full')
    L = opts.reg_pars; % Vary four reg. parameters
    n = length(L);
    [L{:}] = ndgrid(L{end:-1:1});
    L = cat(n+1,L{:});
    RP = fliplr(reshape(L,[],n));
elseif strcmp(opts.reg_vary, 'none')
    RP = [opts.reg_pars0{:}];
end


end