function [] = slider_plot(msh)


S.fh = figure('units','normalized', 'position',[.02,.02,.9,.9],...
    'name','slider_plot');
hr = 1;

xc = [.05, .35, .65];
yc = [.55, .1];
w = .28;
h = .4;
pos = {[xc(1), yc(1),w,h], [xc(2), yc(1),w,h], [xc(3), yc(1),w,h];...
    [xc(1), yc(2),w,h], [xc(2), yc(2),w,h], [xc(3), yc(2),w,h]} ;

velmap = brewermap(20, 'RdBu');
salmap = brewermap(15, 'YlOrBr'); % Brewermap.m (available through Github, use the same code); map = brewermap(N,scheme)
salt = 1;
titles = {'NM', 'OM', 'RWW'};

dir = [-1,1,-1];
S.p = cell(2,3,3);
rot = 'cs';
qs=20;
for vs=1:2
    for ct=1:3
        S.ax{vs,ct} = axes('parent',S.fh,'units', 'normalized', 'position', pos{vs,ct}); % [left bottom width height]
        title([titles{ct},': hour = ', num2str(hr)])
        %title()
        if vs==1
            hold on
            caxis([-1.5,1.5]);
            %disp(S.ax{vs,ct})
            colormap(S.ax{vs,ct}, flipud([103,0,31; % set a colormap (based on colorbrewer)
                178,24,43;
                214,96,77;
                244,165,130;
                253,219,199;
                247,247,247;
                247,247,247;
                209,229,240;
                146,197,222])/255);
            vel = dir(ct)* msh{ct}.(rot).pars(msh{ct}.p.progfgood(:,msh{ct}.adcp_good(hr)))';
            %vel = dir(ct)* msh{ct}.(rot).pars(msh{ct}.p.progfgood_vec(:,hr,1));
            S.p{vs,ct,1} = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),vel, 'LineStyle', '-', 'LineWidth',.2, 'FaceColor', 'flat');
            shading flat
            S.p{vs,ct,3} = quiver(msh{ct}.N, msh{ct}.Z(:,:,1),...
                msh{ct}.(rot).pars(:,:,hr,2)*qs,...
                msh{ct}.(rot).pars(:,:,hr,3)*qs, 'color', 'k','AutoScale','on');
            S.p{vs,ct,2} = plot((msh{ct}.p.nbed),...       % Plot channel bed
                msh{ct}.p.zbed,'-k', 'LineWidth', 1.15);
            %S.p{vs,ct,4} = plot(msh{ct}.p.nbed([1 end]),[0 0],'b','linewidth',2);
            colorbar;
            hold off
            axis tight
        end
        if vs==2 && salt
            colormap(S.ax{vs,ct}, salmap);
            sal = msh{ct}.CTD_int.Sal_ex((msh{ct}.p.progfgood(:,hr)))';
            %sal = msh{1, 1}.CTD_int.Sal{hr, 1};
            S.p{vs,ct,1} = patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
                msh{ct}.p.Z(:,:,1), sal, 'LineStyle', 'none','FaceColor', 'flat'); 
            colorbar;     
            caxis([0,30]);
            hold on
            S.p{vs,ct,2} = plot((msh{ct}.p.nbed),...       % Plot channel bed
                msh{ct}.p.zbed,'-k', 'LineWidth', 1.15);
            hold off
            axis tight
        end
    end
end
%update(S);


% Slider for slope parameter:
S.mSlider = uicontrol('style','slider',...
    'unit','normalized',...  % "pix" is not smart, we are not in a hurry
    'position',[xc(2), .05, w, .02],...
    'min',1,'max',13,'value', hr,...
    'sliderstep',[1/12,1/12],...
    'callback', {@SliderCB, msh});
%A cell array in which the first element is a function handle.
%Subsequent elements in the cell array are the arguments to pass to the callback function.


guidata(S.fh, S);  % Store S struct in the figure
end

function SliderCB(mSlider, EventData, msh)
% Callback for both sliders
S = guidata(mSlider);  % Get S struct from the figure
S.hr = round(get(mSlider, 'Value'));  % Either 'm' or 'b'
S= update(S, msh);
guidata(mSlider, S);  % Store modified S in figure
end

function S=update(S, msh)
hr = S.hr;
%disp(hr)
%disp(S)
xc = [.05, .35, .65];
yc = [.55, .1];
w = .28;
h = .4;
pos = {[xc(1), yc(1),w,h], [xc(2), yc(1),w,h], [xc(3), yc(1),w,h];...
    [xc(1), yc(2),w,h], [xc(2), yc(2),w,h], [xc(3), yc(2),w,h]} ;

velmap = brewermap(20, 'RdBu');
salmap = brewermap(15, 'YlOrBr'); % Brewermap.m (available through Github, use the same code); map = brewermap(N,scheme)

titles = {'NM', 'OM', 'RWW'};
qs=20;

dir = [-1,1,-1];
rot = 'sec';
salt=1;
for vs=1:2
    for ct=1:3
        S.ax{vs,ct}.Title.String = [titles{ct},': hour = ', num2str(hr)];
        if vs==1
            %hold on
            caxis([-1.5,1.5]);
            %colormap(S.ax{vs,ct},velmap);
            vel = dir(ct)* msh{ct}.(rot).pars(msh{ct}.p.progfgood(:,msh{ct}.adcp_good(hr)))';
            %vel = dir(ct)* msh{ct}.(rot).pars(msh{ct}.p.progfgood_vec(:,1,hr));
            S.ax{vs,ct}.Children(3).CData = vel; % Children(3) = patch
            
            S.ax{vs,ct}.Children(2).UData = msh{ct}.sec.pars(:,:,hr,2)*qs;
            S.ax{vs,ct}.Children(2).VData = msh{ct}.sec.pars(:,:,hr,3)*qs;

            % Children(2) = quiver
            %S.ax{vs,ct}.Children(1).CData = vel; % Children(1) = bed plot
            %S.p{vs,ct,2} = plot((msh{ct}.p.nbed),...       % Plot channel bed
            %    msh{ct}.p.zbed,'-k', 'LineWidth', 1.15); 
            colorbar;
            %hold off
            axis tight
        end
        if vs==2 && salt
            %colormap(S.ax{vs,ct}, salmap);
            sal = msh{ct}.CTD_int.Sal_ex((msh{ct}.p.progfgood(:,hr)))';
            %sal = msh{1, 1}.CTD_int.Sal{hr, 1};
            S.ax{vs,ct}.Children(2).CData = sal;
            %colorbar;     
            caxis([0,30]);
            axis tight
        end
    end
end
end