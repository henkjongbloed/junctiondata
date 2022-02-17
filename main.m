%% Read, filter, plot ADCP data
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
close all % Close all figures
%clc
RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

ds = 1; %Choose which dataset to analyze

DS = {'OMHA14','NMOMNW15'}; DS = DS{ds}; 
BN = {{' HK'; ' OMS'; ' OMN'}, {' NM'; ' OM'; ' NWW'}}; BN = BN{ds};

%% Import preprocessed semi-raw data (ADCP / CTD / H structs)

[adcp, ctd, h] = import_data(RF, DS);

% ctd = hctd(h,ctd,DS); %Postprocess ctd data

%% Data analysis





wl = @(t) interp1(datenum(h.time), h.level, t, 'linear');
tak = 1;

%tak = 1:2;
plots = 0;
for b = tak
    %plot_all(adcp{b})
    [V{b}, U{b}] = processVel(adcp{b}, h, BN{b}); % V = VMADCP objects, U = generic velocity data (plus more)
    S{b} = processSal(V{b}, U{b}, ctd{b}, DS); %S = Salinity data
    if plots
        try
            plot_all(V{b})
            %plot_geom(U{b})
        catch
            disp('Some error in plotting')
        end
    end
end

% for b=tak
%     F{b} = processFlux(U{b}, S{b});
% end
%%
figure;
plot_geom_all(U, S, h)

%% Study tidal phase and amplitudes
%create preliminary plot to view results



for j = tak
    for i = 1
        subplot(length(tak), 1,sub2ind([length(tak), 1], j,i ))
        U{j}.mesh_mean.plot(cell2mat(U{j}.n_bvels{1}))
        %caxis([min(min(U{j}.tid_pars(:,d))), max(max(U{j}.tid_pars(:,d)))])
        colormap(gca, velmap)
        colorbar;
    end
end
% Tidal

figure;
% suptitle('Flow')
tit = {'Branch 1, Const 0x', 'Branch 2, Const 0x', 'Branch 3, Const 0x';...
    'Branch 1, A2', 'Branch 2, A2x', 'Branch 3, A2x';...
    'Branch 1, Phi2x', 'Branch 2, Phi2x', 'Branch 3, Phi2x';...
    'Branch 1, Const 0y', 'Branch 2y, Const 0', 'Branch 3, Const 0y';...
    'Branch 1, Const 0z', 'Branch 2, Const 0z', 'Branch 3, Const 0z'};
%U{1, 1}.T.velocity_model  
%d = 3;

cnew = [velmap ;flipud(velmap)];
for j = tak
    for i = 1:size(U{j}.tid_pars,2)
        subplot(size(U{j}.tid_pars,2), length(tak), sub2ind([length(tak), size(U{j}.tid_pars,2)], j,i))
        U{j}.mesh_mean.plot(U{j}.tid_pars(:,i))
        title(tit{i,j})
        if i==3
        %caxis([min(min(U{j}.tid_pars(:,d))), max(max(U{j}.tid_pars(:,d)))])
            colormap(gca, cnew)
            caxis([-pi, pi])
        else
            colormap(gca, velmap)
        end
        colorbar;
    end
end

%% 
%create preliminary plot to view results
figure;
% suptitle('Flow')    
d = 2;

ct = 1;

for j = tak
    for i = 1:numel(U{j}.mesh)
        subplot(length(tak),numel(U{j}.mesh),sub2ind([length(tak),13], j,i))
        U{j}.mesh(i).plot(U{j}.vel(:,d))
        caxis([min(min(U{j}.vel{i}(:,d))), max(max(U{j}.vel{i}(:,d)))])
        colormap(gca, velmap)
        colorbar;
    end
end

%%
figure;
suptitle('Salt')
for i =1:12
    for j = tak
        ind = sub2ind([length(tak),13], j,i);
        subplot(13,length(tak),ind)
        U{j}.mesh(i).plot(S{j}.S{i}(:,1))
        colormap(gca, salmap)
        colorbar;
    end
end




%% Salt flux plots
XTL = {'F_{00}', 'F_{01}^*', 'F_{02}^*', 'F_{03}^*', 'F_{04}^*',...
    'F_{10}', 'F_{11}^*',...
    'F_{20}', 'F_{21}^*',...
    'F_{30}', 'F_{31}^*', 'F_{32}^*', 'F_{33}^*', 'F_{34}^*', 'F_{35}', 'F_{36}^*',...
    'F_{40}', 'F_{41}', 'F_{42}'};


dims = {' us', ' vs'};

figure;
dim = 1:2;
for b = tak
    for d = dim
        subplot(length(dim),length(tak),sub2ind([length(tak),length(dim)], b,d))
        hold on
        for i = 1:5
            plot(U{1, b}.mesh_mean.n_middle, F{1,b}.Loc{i,d})
            grid on
            legend('Subtidal Adv.','Tidal-Sal. CC','Stokes drift','Tidal Sloshing','Vert. Shear disp.')
        end
        title(strcat('Kjerfve Decomposition', dims{d}, ' at ', BN{b}))
    end
end
        
figure;
dim = 1:2;
for b = tak
    for d = dim
        subplot(length(dim), length(tak), sub2ind([length(tak), length(dim)], b, d))
        bar(F{1, b}.Lat(:,d))
        title(strcat('Lateral Decomposition', dims{d}, ' at ', BN{b}))
        xticks(1:19)
        xticklabels(XTL)
        grid on
    end
end

%plotKjerfve(TA, mesh_mean)

%plot19(TA)




