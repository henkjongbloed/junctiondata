function [row, col, val] = dom2rowcol(par_names_tot, doma, c, adj, const, row, val)

if strcmp(const, 'M0A')
    if doma == 0
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',c,': d^1u/dx^1: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(1),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(3),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': w0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': w0: M0A')))];
    elseif doma ==1
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',c,': d^1u/dx^1: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(3),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': w0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': w0: M0A')))];
        row(4) = []; val(4) = [];
    elseif doma ==2
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',c,': d^1u/dx^1: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(3),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': w0: M0A')))];
        row([2,4,6,8]) = []; val([2,4,6,8]) = [];
    elseif doma ==3
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',c,': d^1u/dx^1: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(1),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(3),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': w0: M0A')))];
        row([2,6,8]) = []; val([2,6,8]) = [];
    elseif doma ==4
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',c,': d^1u/dx^1: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(1),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': w0: M0A')))];
        row([2,5,6,8]) = []; val([2,5,6,8]) = [];
    elseif doma ==5
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',c,': d^1u/dx^1: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(1),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': w0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(4),': w0: M0A')))];
        row(5) = []; val(5) = [];
    elseif doma ==6
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',c,': d^1u/dx^1: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(1),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': w0: M0A')))];

        row([3,5,7,9]) = []; val([3,5,7,9]) = [];
    elseif doma ==7
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',c,': d^1u/dx^1: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(1),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(3),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': w0: M0A')))];
        row([3,7,9]) = []; val([3,7,9]) = [];
    else
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',c,': d^1u/dx^1: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': u0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(3),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': v0: M0A'))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',adj(2),': w0: M0A')))];
        row([3,4,7,9]) = []; val([3,4,7,9]) = [];
    end

else % Tidal constituent introduces two extra terms in the continuity eq.

    if doma == 0
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dx^1: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': w0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': w0: ',const)))];
    elseif doma ==1
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dx^1: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': w0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': w0: ',const)))];
        row(4) = []; val(4) = [];
    elseif doma ==2
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dx^1: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': w0: ',const)))];
        row([2,4,6,8]) = []; val([2,4,6,8]) = [];
    elseif doma ==3
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dx^1: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': w0: ',const)))];
        row([2,6,8]) = []; val([2,6,8]) = [];
    elseif doma ==4
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dx^1: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': w0: ',const)))];
        row([2,5,6,8]) = []; val([2,5,6,8]) = [];
    elseif doma ==5
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dx^1: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': w0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(4),': w0: ',const)))];
        row(5) = []; val(5) = [];
    elseif doma ==6
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dx^1: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': w0: ',const)))];

        row([3,5,7,9]) = []; val([3,5,7,9]) = [];
    elseif doma ==7
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dx^1: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(1),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': w0: ',const)))];
        row([3,7,9]) = []; val([3,7,9]) = [];
    else
        col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',c,': d^1u/dx^1: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': u0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(3),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': v0: ',const))) , ...
            find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',adj(2),': w0: ',const)))];
        row([3,4,7,9]) = []; val([3,4,7,9]) = [];
    end
    col = [col, find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',c,': d^1u/dx^1: M0A'))), ...
        find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',c,': d^1v/dy^1: M0A')))]; % Two extra terms
end
end