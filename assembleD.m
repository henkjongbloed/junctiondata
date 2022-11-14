function [D, IM] = assembleD(LBS)

par_names = LBS.velocity_model.names;
Np = sum(LBS.velocity_model.npars);
NNp = Np*LBS.mesh.ncells;
idx = 1;
w = ones([Np,1]);
IM = zeros(LBS.mesh.ncells);
for j = 1:LBS.mesh.ncells % loop trough every cell
    for i=1:Np
        par_names_tot{1,idx} = sprintf('cell %i: %s',j, par_names{1,i});       
        idx = idx + 1;
    end
    [neighbors(j,:), dom(j)] = LBS.mesh.get_neighbors(j);
    nbreal = neighbors(j,~isnan(neighbors(j,:)));
    nbreal = nbreal(nbreal>j);
    IM(j,nbreal) = 1;
end

IM = IM + IM';

% Apply enhanced regularization for small singular value features

for i = 1:Np
%     if contains(par_names{i}, 'M4')
%         w(i) = w(i)*10;
%     end
    if contains(par_names{i}, 'v') || contains(par_names{i}, 'w')
        w(i) = w(i)*10;
    end
    if contains(par_names{i}, 'dy') || contains(par_names{i}, 'dx')
        w(i) = w(i)*100;
    end
end
W = sparse(diag(repmat(w,LBS.mesh.ncells,1)));
D1 = speye(NNp);
rows = []; cols = []; vals = [];
for c = 1:LBS.mesh.ncells %rows
    adj = neighbors(c,:);
    adj = adj(~isnan(adj)); nnb = length(adj);
    for nb = 1:nnb
        col = ((adj(nb)-1)*Np+1):(adj(nb)*Np);
        row = (c-1)*Np+1:c*Np;
        val = -1/nnb*ones(1,Np);    
        rows = [rows row];
        cols = [cols col];
        vals = [vals val];
    end
%     [row, col, val] = dom2rowcol(par_names_tot, doma, c, adj, 'M0A', row, val);

%     row_idx = row_idx + 1;
end
D = D1 + sparse(rows, cols, vals, NNp, NNp);

D = W*D;

end