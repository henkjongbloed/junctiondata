function addPaths(RF, DS, varargin)

% SET PATHS: [SET YOUR OWN]
%RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\';
for ca=1:numel(varargin)
    if strcmp(varargin{ca}, 'old')
        rmpath(strcat(RF,'Tools\adcptools'));
        addpath(strcat(RF,'Tools\adcptools_r94'));
    elseif strcmp(varargin{ca}, 'new')
        rmpath(strcat(RF,'Tools\adcptools_new'));
        addpath(strcat(RF,'Tools\adcptoolsGit'));
    end
end
addpath(strcat(RF,'Tools\loess-master'));
addpath(strcat(RF,'Tools\T_Tide'))
%addpath(strcat(RF,'Tools\example_adcp_processing\muaramuntai25Aug'));
addpath(strcat(RF,'Tools\BrewerMap-master'));
addpath(strcat(RF,'Tools\subaxis'));
addpath(strcat(RF,'Data\processedData'));

%addpath(strcat(RF,'processedData'));

addpath(strcat(RF,'Data\',DS));
addpath(strcat(RF,'Data\MerelFiles\_processed'));

%addpath(strcat(RF,'DataAnalysis\ProcessedData'));

%addpath(strcat(RF,'Tools\BrewerMap-master'));

cd(strcat(RF,'DataAnalysis'))

end
