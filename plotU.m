function plotU(TA, mesh, var)

br = length(TA);
velmap = brewermap(20, 'RdBu');

for b = 1:br
    if strcmp(var,'u')
        u = TA{b}.u;
    else
        u = TA{b}.s;
    end
    nc = mesh{b}.ncells ;
    siz = size(mesh{b}.col_to_mat);
    nrows = siz(1); 
    ncols = siz(2);
    for c = 1:5
        subplot(5,br,sub2ind([br,5],b, c))
        %disp(sub2ind([5,br],c,b))            
        if c == 1
            utemp = u{1,1}*ones(nc,1);
            mesh{b}.plot(utemp)

        elseif c==2
            utemp = repmat(u{2,1},nrows,1);
            mesh{b}.plot(utemp(mesh{b}.mat_to_cell))
        elseif c==3
            plot(u{3,1})
        elseif c==4
            contourf(u{4,1})
        else
            utemp = mean(u{5,1},1,'omitnan');
            mesh{b}.plot(utemp(mesh{b}.mat_to_cell))
        end
        if c==1 || c==2 || c==5
            shading flat
            %axis equal
            set(gca,'ylim', [min(mesh{b}.zb_all), mesh{b}.water_level])
            set(gca,'xlim', [mesh{b}.nw(1), mesh{b}.nw(2)] )
        elseif c==3 || c==4
        end
            
        %caxis([-1.5,1.5])
        colormap(gca, velmap)
        colorbar;
        title(strcat(var, ': Component =  ',num2str(c)));

    end
end