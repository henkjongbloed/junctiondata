function u = p2u(p, names, eta, mesh, t)

% Convert tidal parameter vector to velocities on a mesh: Evaluation for a
% sequence of temporal points t

% Now: Only M2 and M4 compatible

% Output u will be array of size Np/(1 + 2*Nc) x length(t)

nc = 2*length(eta.omega)+1;
sizu = size(p,1)./nc;
u = zeros(sizu, length(t));

ind = {find(contains(names, 'M0A')),...
        find(contains(names, 'M2A')),...
        find(contains(names, 'M2B')),...
        find(contains(names, 'M4A')),...
        find(contains(names, 'M4B'))};

% Because of ordering of parameter vector, we do not need to explicitly
% find indices
tidvec = {ones(size(t)), cos(eta.omega(1)*t), sin(eta.omega(1)*t),...
     cos(eta.omega(2)*t), sin(eta.omega(2)*t)};

for nam = 1:sizu
    locp = nc*(nam-1)+1;
    u(nam, :) = p(locp,1).*tidvec{1};
    for c = 1:nc
        u(nam, :) = u(nam, :) + p(locp + c, 1)*tidvec{c} + p(locp + c, 1)*tidvec{c};
    end
end
    






% for i = 1:length(ind)
%     u(ind{i}) = p*;
% end


end