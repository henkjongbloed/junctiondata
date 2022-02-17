function [x,y,z] = cells2mat(U, opt)
%Process velocities

if strcmp(opt, 'vel')
    cm = cell2mat(U.vel');
    
    x = cm(:,1:3:end);
    y = cm(:,2:3:end);
    z = cm(:,3:3:end);
    
elseif strcmp(opt,'sal')
    x = cell2mat(U.S);
    y = x;
    z = x;
end

end