function [row, col, val] = dom2rowcol2(par_names_tot, doma, c, adj, const, dn, dsig)
% Helper function for the extended consistency matrix assembly
% (assembleC4.m)

% First y derivatives, then sigma derivatives in order of u, v, w
row = [0,0,0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5]; % Row incremental index, always the same
val = [1, -1/(dn), 1/(dn),...
    1, -1/(dn), 1/(dn),...
    1, -1/(dn), 1/(dn),...
    1, -1/(dsig), 1/(dsig),...
    1, -1/(dsig), 1/(dsig),...
    1, -1/(dsig), 1/(dsig)]; % Value, always the same (but look at order of magnitudes difference between the values...)

if doma == 0
    col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': d^1u/dy^1: ', const)))];
elseif doma == 1
    col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': d^1u/dy^1: ', const)))];
elseif doma == 2
    col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': d^1u/dy^1: ', const)))];
elseif doma == 3
    col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': d^1u/dy^1: ', const)))];
elseif doma == 4
    col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': d^1u/dy^1: ', const)))];
elseif doma == 5
   col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': d^1u/dy^1: ', const)))];
elseif doma == 6
   col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const)))];
elseif doma == 7
    col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const)))];
elseif doma == 8
    col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': w0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': u0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1v/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': v0: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1w/dsig^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': d^1u/dy^1: ', const))) , ...
        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dy^1: ', const)))];
end

end