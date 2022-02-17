function plotKjerfve(TA, mesh)
BN = {' NM', ' OM', ' NWW'};
br = length(TA);
cr = 5;
figure;
for b = 1:br
    for c = 1:cr
        subplot(cr,br,sub2ind([br,cr],b, c))
        if c ==1 || c==2
            hold on
            plot(mesh{b}.n_middle, TA{b}.F{1,c})
            plot(mesh{b}.n_middle, TA{b}.F{2,c})
            plot(mesh{b}.n_middle, TA{b}.F{3,c})
            plot(mesh{b}.n_middle, TA{b}.F{4,c})
            plot(mesh{b}.n_middle, TA{b}.F{5,c})
            %plot(mesh_mean.n_middle, TA.Fx)
            %plot(mesh_mean.n_middle, TA.FD)
        elseif c==3
            plot(repmat(mesh{b}.n_middle,13, 1) + 20*TA{b}.vDA, TA{b}.uDA)
            title(strcat('Tidal Ellipses (raw) at ', BN{b}))
        elseif c==4
            hold on
            plot(mesh{b}.n_middle,TA{b}.A{1,1}(1,:))
            plot(mesh{b}.n_middle,TA{b}.A{1,1}(2,:))
            plot(mesh{b}.n_middle,TA{b}.A{1,2}(1,:))
            plot(mesh{b}.n_middle,TA{b}.A{1,2}(2,:))
            legend('AM2', 'AM4', 'PhiM2', 'PhiM4')
            title('Flow: Tidal Fit')
        elseif c==5
            hold on
            plot(mesh{b}.n_middle,TA{b}.A{3,1}(1,:))
            plot(mesh{b}.n_middle,TA{b}.A{3,1}(2,:))
            plot(mesh{b}.n_middle,TA{b}.A{3,2}(1,:))
            plot(mesh{b}.n_middle,TA{b}.A{3,2}(2,:))
            legend('AM2', 'AM4', 'PhiM2', 'PhiM4')
            title('Salt: Tidal Fit')
            if c==1
                title(strcat('Kjerfve Decomposition at: ', BN{b}))
                if b==br
                    legend('F0: Adv', 'F1: Tide/Sal', 'F2: Stokes', 'F3: Sloshing', 'F4: Shear')
                end
            end
        end
    end
end
