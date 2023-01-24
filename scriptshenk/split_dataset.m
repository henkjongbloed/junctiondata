function training_idx = split_dataset(opts)

training_idx = ones(size(opts.cell_idx));

if strcmp(opts.cv_mode, 'none')
    training_idx = ones(size(opts.cell_idx));
elseif strcmp(opts.cv_mode, 'random')
    tp = opts.training_perc;
    rand0 = rand(size(training_idx));
    training_idx = (rand0 <= tp);
elseif strcmp(opts.cv_mode, 'omit_cells')
    oc = opts.omit_cells;
    for occ = 1:length(oc)
        training_idx(opts.cell_idx==oc(occ)) = 0;
    end
elseif strcmp(opts.cv_mode, 'omit_time') % to be implemented
end

end