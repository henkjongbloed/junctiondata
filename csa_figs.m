function csa_figs(didx, bidx, D, flow, AF, tak_name)
% Selection of data
c = {'k', 'k', 'k'};
ls = {'-', '--', '-.'};

% AF{di}{bi}{vi}{ri, c}
tD = D{didx(1)}{bidx(1)}.X.t  ; % days
for di = didx 
    for bi = bidx
        t0 = flow{di}{bi}.solver.adcp.time(1); % start of experiment  dataset di branch bi
        trealn = datenum(t0) + tD;%datenum unit = 1 -> datetime diff = 1 day.
        treal{di}{bi} = datetime(trealn, 'ConvertFrom', 'datenum');
    end
end

for di = didx    
    %[fig, TL] = prep_fig([20,9], [3,1]);

    %waterlevel - same for all bracnhes
%     opts.ylim = "sym";
opts.ylim = 0;
    ax{1} = prep_ax(); 

    axis tight
    ax{1}.XTickLabel = [];
    ylabel("m", "FontName","BookAntiqua", "Interpreter","latex")
    yyaxis left
    hold on
    %treal_c{di} = sort([treal{di}{:}]);

%     yyaxis right
    for bi = bidx
%         treal_c{di} = [treal{di}{:}]; %flow{1, 1}{1, 1}.solver.adcp.water_level_object.model
        eta{di}{bi} = flow{di}{bi}.solver.adcp.water_level_object.get_water_level(treal{di}{bi});
        eta_model{di}{bi} = flow{di}{bi}.solver.adcp.water_level_object.get_water_level_model(treal{di}{bi});
        plot(treal{di}{bi}, eta{di}{bi}, c{1}, 'LineWidth',1)
        plot(treal{di}{bi}, eta_model{di}{bi}, c{1}, 'LineWidth',1)
%         plot(treal{di}{bi}, D{di}{bi}.A)
    end
    hold off
%     plot(treal, eta{1,di}, c{1}, 'LineWidth',1)
    title("Waterlevel", "FontName","BookAntiqua")
    fix_ax(ax{1}, opts)

    %flow
    ax{2}= prep_ax();hold on
    for bi = bidx
        plot(treal{di}{bi} , AF{di}{bi}{1}{1}{2}, 'Color', c{2}, 'LineStyle', ls{bi}, 'LineWidth',1)
    end
    hold off
    ax{2}.XTickLabel = [];
    ylabel("m/s", "FontName","BookAntiqua", "Interpreter","latex")
    legend(tak_name{di}{:})
    axis tight
    title("Flow velocity", "FontName","BookAntiqua")
    %sumu = D{1,di}.A.* AF{1,di}{1,2} + D{2,di}.A.* AF{2,di}{1,2}- D{3,di}.A.* AF{3,di}{1,2};
    %plot(taxis , sumu, 'Color', c{2}, 'LineStyle', ls{bi})


    %salt
    fix_ax(ax{2}, opts)
    ax{3}  = prep_ax();hold on
    for bi = bidx
        plot(treal{di}{bi}, AF{di}{bi}{2}{1}{2}, 'Color', c{3}, 'LineStyle', ls{bi}, 'LineWidth',1)

    end
    hold off
    ylabel("psu", "FontName", "BookAntiqua", "Interpreter", "latex")
    axis tight
    title("Salinity", "FontName","BookAntiqua")
    legend(tak_name{di}{:})
    opts.ylim = "pos";
    ax{3} = fix_ax(ax{3}, opts);
    TL.TileSpacing = 'compact';
    TL.Padding = 'compact';
    %linkaxes([ax{:}],'x')
    TL.XLabel.String = 't [h]';
%     t.YLabel.String = 'z';
    TL.XLabel.Interpreter = 'latex';
    TL.YLabel.Interpreter = 'latex';
    
end
end