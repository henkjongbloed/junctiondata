%slider_main

S.x   = 0:.01:1;  % For plotting.         
S.fcn = @(x, m, b) m * x + b; % in case of msh, only 13 possible input data (hours)
S.m   = 1;
S.b   = 0;

slider_plot(S)