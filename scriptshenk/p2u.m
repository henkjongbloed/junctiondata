function u = p2u(p, names, eta, mesh, t)

% Convert tidal parameter vector to velocities on a mesh: Evaluation for a
% sequence of temporal points t

% Now: Only M2 and M4 compatible

% Output u will be array of size Np/(1 + 2*Nc) x length(t)


u = zeros(size(p,1), length(t));
ind = {find(contains(names, 'M0A')),...
        find(contains(names, 'M2A')),...
        find(contains(names, 'M2B')),...
        find(contains(names, 'M4A')),...
        find(contains(names, 'M4B'))};

% Because of ordering of parameter vector, 
%for nam = 1:length(names)
    

tidfun = {1, @(t) cos(eta.omega(1)*t), @(t) sin(eta.omega(1)*t),...
    @(t) cos(eta.omega(2)*t), @(t) sin(eta.omega(2)*t)};





% for i = 1:length(ind)
%     u(ind{i}) = p*;
% end


end