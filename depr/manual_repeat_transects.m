function ef_out = manual_repeat_transects(V)

fi = V.fileid;

nr = max(fi); %assuming that fi starts from 1 to nr
ef(nr,1) = EnsembleFilter; % initialize output ensemble filter array
for rt = 1:nr
    good = (fi==rt*ones(size(fi)));
    ef(rt) = EnsembleFilter(~good);
end

ef_out = ef;

end
    