function zb = n2zb(n, m)
% Get bed elevation from n coordinate, mesh m.

zb = interp1(m.nb_all, m.zb_all, n, 'linear');
end