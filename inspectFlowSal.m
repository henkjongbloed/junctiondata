function inspectFlowSal(msh)
% Create figure window and components

fig = uifigure('Name','Inspect Flow and Salinity');
% ct = 1;
% hr = 1;
% dir = [-1,1,-1];
% velmap = brewermap(20, 'RdBu');
% 
% for ct = 1:3
%     titles = {'NM', 'OM', 'RWW'};
%     sp(ct) = subplot(3, 1, ct);
%     title([titles{ct},': hour = ', num2str(hr)])
%     hold on
%     caxis([-1.5,1.5]);
%     vel = dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood(:,msh{ct}.adcp_good(hr)))';
%     p(ct) = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),vel, 'LineStyle', '-', 'LineWidth',.2, 'FaceColor', 'flat');
%     q(ct) = quiver(msh{ct}.N, msh{ct}.Z(:,:,1),...
%         msh{ct}.sec.pars(:,:,hr,2),...
%         msh{ct}.sec.pars(:,:,hr,3), 'color', 'k', 'LineWidth',1.5,'ShowArrowHead','on','AutoScale','on');
%     b(ct) = plot((msh{ct}.p.nbed),...       % Plot channel bed
%         msh{ct}.p.zbed,'-k', 'LineWidth', 1.15); colormap(velmap); 
%     cb = colorbar;
%     hold off
%     axis tight
% end
i = 0;
p = plot(1:10, i.*[1:10])

sld = uislider(fig,'ValueChangedFcn',@(sld,event) updateHour(sld));

end

% Create ValueChangedFcn callback
function hr = updateHour(sld)

p.YData = sld.*[1:10];
end