function plotPe(dat)

rp = dat.opts.reg_pars;
for i = 1:size(dat.Pe, 1)
    Pe(:,:,i) = reshape(dat.Pe(i,:)',length(rp{3}), length(rp{1}),[]);
end
% regP = dat.regP;
if strcmp(dat.opts.reg_vary, 'coupled')
    [X,Y] = meshgrid(symlog(rp{1}), symlog(rp{3}));
%     X = unique(X','rows')';
%     X = unique(X);
%     Y = unique(Y,'rows')';
end
tit = {'Total error', 'Training Error', 'Test error', 'Continuity: Inside cell', 'Continuity: Inter-cell', 'Coherence', 'Consistency'};
for fi = 1:length(tit)
    subplot(3,3,fi)
    contourf(X,Y,Pe(:,:,fi)./Pe(1,1,fi) ,20,'LineWidth',1)
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    xlabel('lc')
    ylabel('ld')
    title(tit{fi})
    colormap("jet")
% 
%     xlabel('sB')
%     ylabel('l_c (log)')
%     title(sprintf('||p - p_e||. l_p = %1.1d', Lp(lp_ind(fi))))
    colorbar;
%     colormap(flipud(hot));     caxis(cap)
end

%     subplot(2,4,fi+4)
%     contourf(sB, Lc, log10(squeeze(solB(:,lp_ind(fi), :))),20,'LineWidth',1)
%     xlabel('sB')
%     ylabel('lc (symlog)')
%     title(sprintf('%s%d', 'MSE B estimation. lp:', Lp(lp_ind(fi))))
%     colorbar;
%     colormap(flipud(hot));     caxis(cab)

end