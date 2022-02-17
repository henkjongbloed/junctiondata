%compare msh data

clearvars
close all
clc
% SET PATHS: [SET YOUR OWN]
Path_ADCPtools      = 'C:\Users\jongb013\OneDrive - WageningenUR\PhDHenkJongbloed\2. Programming\DataTools\Tools\ADCPTools30-06-20';
Path_BrewerMap      = 'C:\Users\jongb013\OneDrive - WageningenUR\PhDHenkJongbloed\2. Programming\DataTools\Tools\BrewerMap-master';
Path_MATLAB         = pwd;%'C:\Users\Judith\Documents\Stage_HbR\_Scripts\MATLAB';

addpath('C:\Users\jongb013\OneDrive - WageningenUR\PhDHenkJongbloed\2. Programming\DataTools\DataHartelkanaal\14sep2015zout')
addpath(Path_ADCPtools)
addpath(Path_MATLAB)
addpath(Path_BrewerMap)

mshSalt = load('MSH_20150914_v1.mat','msh');
mshVelo = load('WURData_processed.mat','msh');

mshSalt = mshSalt.msh;
mshVelo = mshVelo.msh;

% Conclusion: mesh and thus all other quantities differ considerably. 