function F = Kjerfve(h, u, v, s)

h1 = h{3};
h0 = h{1} + h{2};

hs = repmat(h0, length(h1), 1) + h1;

%hs = h{1} + h{2} + h{3};

shearTermx = mean(u{5,1}.*s{5,1},3,'omitnan');
shearTermy = mean(v{5,1}.*s{5,1},3,'omitnan');
ny = numel(u{1,2});

F{1,1} = h0.*u{1,2}.*s{1,2};
F{2,1} = u{1,2}.*mean(repmat(h1,1,ny).*s{2,2},1,'omitnan');
F{3,1} = s{1,2}.*mean(repmat(h1,1,ny).*u{2,2},1,'omitnan');
F{4,1} = mean(hs.*u{2,2}.*s{2,2}, 1, 'omitnan');
F{5,1} = mean(hs.*shearTermx);

F{1,2} = h0.*v{1,2}.*s{1,2};
F{2,2} = v{1,2}.*mean(repmat(h1,1,ny).*s{2,2},1,'omitnan');
F{3,2} = s{1,2}.*mean(repmat(h1,1,ny).*v{2,2},1,'omitnan');
F{4,2} = mean(hs.*v{2,2}.*s{2,2}, 1, 'omitnan');
F{5,2} = mean(hs.*shearTermy);



end