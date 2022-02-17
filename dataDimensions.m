
% transect = 1; %1, 2 or 3
% In our data, this is the case:
% CS = 3;% - number of cross-section
% E = [length(adcp201408{1}.FileNumber), length(adcp201408{2}.FileNumber), length(adcp201408{3}.FileNumber)];% - number of ensembles
% R = [size(msh201408{1}.pars,1), size(msh201408{2}.pars,1), size(msh201408{3}.pars,1)];%- maximum number of rows in output mesh
% C = [size(msh201408{1}.pars,2), size(msh201408{2}.pars,2), size(msh201408{3}.pars,2)];%- number of columns in output mesh
% T = [size(msh201408{1}.pars,3), size(msh201408{2}.pars,3), size(msh201408{3}.pars,3)];%- number of time steps in output
% P = [size(msh201408{1}.pars,4), size(msh201408{2}.pars,4), size(msh201408{3}.pars,4)];%- number of fitted parameter (i.e. sum of the parameters for
%                                 %x,y and z components of the velocity)
% Q = [size(msh201408{1}.p.X,2), size(msh201408{2}.p.X,2), size(msh201408{3}.p.X,2)];%- number of cells in the cross-section