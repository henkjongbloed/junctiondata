function [lf, ls, ures, sres, evres, dname, bname, flow, salt, F, D, DF, AF, NF, SF, adcp, ctd, h] = preallo(l)
% l = 10;
addpath("helpers\", "saved\", "plot\")
warning("off","MATLAB:subscripting:noSubscriptsSpecified")

lf = l*[10,10,1,1,10];                  %Reg pars flow
ls0 = 1;
ls = l*ls0*[1,1];                           %Reg pars salt

ures = [50, 21];                          %Mesh resolution flow
sres = [10, 21];                          %Mesh resolution salinity

evres = [26, ures(1), ures(2)];           %Time x Lateral x Vertical - decomp resolution.

dname = {'OMHA14','NMOMNW15'};
bname = {{' Hartel Canal'; ' Old Meuse South'; ' Old Meuse North'},...
    {' New Meuse'; ' Old Meuse'; ' New Waterway'}};

%% Get solution
flow = cell(2,1);
salt = cell(2,1);
F = cell(2,1); % 3 variables: flow u, salt s, flux us
D = cell(2,1); % indep of reg

DF = cell(2,1);
AF = cell(2,1);

NF = cell(2,1);
SF = cell(2,1);

adcp = cell(2,1);

ctd = cell(2,1);
h = cell(2,1);

end