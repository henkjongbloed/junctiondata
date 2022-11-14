function C = assembleC(LBS,idx)

% Function that assembles cell-based continuity equation
cmesh = LBS.mesh;
ZB = LBS.bathy;
T = LBS.velocity_model;
eta = LBS.adcp.water_level_object;
B = 1;
L = 1;
sig = cmesh.sig_center(idx);
D0 = eta.parameters(1) - LBS.mesh.zb_middle(LBS.mesh.col_to_cell(idx));
D0x = 0; % Assumption of zero gradient in x-direction
dy = .1;
D0y = 0;%-1/(2*dy)*(cmesh.zb.zbI(cmesh.r0sig(1,cell_idx) + dy) - cmesh.zb.zbI(cmesh.r0sig(1,cell_idx) - dy));
% We consider B = L = 1 for now (no scaling of horizontal coordinate)
C = zeros([1+2*numel(eta.constituents), sum(T.npars)]);
% ind = zeros([3, numel(names)]);
% First: subtidal equation

continuity_possible = 0;
if (T.s_order(1) > 0) && (T.n_order(2) > 0) && (T.sigma_order(3) > 0 || T.z_order(3) > 0)
    continuity_possible = 1;
end
if continuity_possible
    ind0 = [find(strcmp(T.names, 'd^1u/dx^1: M0A')) , ...
        find(strcmp(T.names, 'd^1u/dsig^1: M0A')) , ...
        find(strcmp(T.names, 'd^1v/dy^1: M0A')) , ...
        find(strcmp(T.names, 'd^1v/dsig^1: M0A')) , ...
        find(strcmp(T.names, 'd^1w/dsig^1: M0A'))];

    terms0 = [B*D0, B*(1-sig)*D0x, L*D0, L*(1-sig)*D0y, L*B];

    C(1,ind0) = terms0;


    % Second: Tidal equations for each constituent
    for i=1:numel(T.constituentsU) %Fundamental choice: All constituents are the same!!
        const = T.constituentsU{i};

        ind{i,1} = [find(strcmp(T.names, sprintf('%s%s%s', 'd^1u/dx^1: ', const,'A'))) , ...
            find(strcmp(T.names, 'd^1u/dx^1: M0A')) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1u/dsig^1: ', const,'A'))) ,...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1v/dy^1: ', const,'A'))) , ...
            find(strcmp(T.names, 'd^1v/dy^1: M0A')) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1v/dsig^1: ', const,'A'))) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1w/dsig^1: ', const,'A')))];

        terms{i,1} = [B*D0, B*eta.parameters(2*i), B*(1-sig)*D0x, L*D0, B*eta.parameters(2*i), L*(1-sig)*D0y, L*B];
        C(2*i,ind{i,1}) = terms{i,1};


        ind{i,2} = [find(strcmp(T.names, sprintf('%s%s%s', 'd^1u/dx^1: ', const,'B'))) , ...
            find(strcmp(T.names, 'd^1u/dx^1: M0A')) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1u/dsig^1: ', const,'B'))) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1v/dy^1: ', const,'B'))) , ...
            find(strcmp(T.names, 'd^1v/dy^1: M0A')) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1v/dsig^1: ', const,'B'))) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1w/dsig^1: ', const,'B')))];

        terms{i,2} = [B*D0, B*eta.parameters(2*i+1), B*(1-sig)*D0x, L*D0, B*eta.parameters(2*i+1), L*(1-sig)*D0y, L*B];
        C(2*i+1,ind{i,2}) = terms{i,2};
    end
end
end
