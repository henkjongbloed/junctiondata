% Fri Oct 17 11:09:36 CEST 2014
% Karl Kastner, Berlin

% Adapted by Judith Poelman, 9-9-2019
function  addPaths_old(ROOTFOLDER)
if ischar(ROOTFOLDER)
% userpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src',filesep]);
addpath([ROOTFOLDER, filesep, 'Scripts\adcptools-code-r80-trunk']);
% addpath([ROOTFOLDER, filesep, '_Scripts\MATLAB_Kastner\example_Karl']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src\lib']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/backscatter']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/instrumentation']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/instrumentation/adcp']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/interpolation']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/linear-algebra']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/linear-algebra/coordinate-transformation']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/mathematics']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/mathematics/numerical-methods']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/mathematics/numerical-methods/finite-difference']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/mathematics/numerical-methods/integration']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/mathematics/signal-processing']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/mathematics/statistic']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/mathematics/statistic/moment-statistics']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/mathematics/statistic/order-statistics']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/mathematics/statistic/resampling-statistics']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/matlab']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/plot']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/regression']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/strings']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/system']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib/velocity-profile']);
% addpath([ROOTFOLDER,filesep,'_Scripts\MATLAB_Kastner\src/lib\mathematics\geometry']);
% %addpath([ROOTFOLDER, filesep,'_Scripts\20140812_Botlek']); 
% addpath([ROOTFOLDER, filesep, '_Scripts\MATLAB_Kastner\src/master1']);
% addpath([ROOTFOLDER, filesep, '_Scripts\MATLAB_Kastner\src/master1/fem']);
% addpath([ROOTFOLDER, filesep, '_Scripts\Discharge']);
addpath([ROOTFOLDER, filesep, '\Data\master']); 
addpath([ROOTFOLDER, filesep, '\Data\slave']);  
addpath([ROOTFOLDER, filesep, '\Data']); 
%addpath([ROOTFOLDER, filesep, '_Scripts/MATLAB_HbR']); 
addpath([ROOTFOLDER, filesep, 'Scripts\altmany-export_fig-dd9397c']);   


set(0,'DefaultAxesFontName', 'CMU Serif');
set(0,'DefaultTextFontName', 'CMU Serif')
end
end


