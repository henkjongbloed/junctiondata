function plotSe(dat, noise_ind)

rp = dat.opts.reg_pars;
SA = reshape(dat.avg_err,length(rp{3}), length(rp{1}),[]);
if strcmp(dat.opts.reg_vary, 'coupled')
    [X,Y] = meshgrid(symlog(rp{1}), symlog(rp{3}));
end
figure;
% tit = {'Total error', 'Training Error', 'Test error', 'Continuity: Inside cell', 'Continuity: Inter-cell', 'Coherence', 'Consistency'};
for fi = 1:size(SA,3)
    subplot(3,3,fi)
    contourf(X,Y, SA(:,:,fi)./SA(1,1,fi) ,20,'LineWidth',1)
    xlabel('lc')
    ylabel('ld')
%     title(tit{fi})
    colormap("jet")
    colorbar;
end
figure;
SAP = dat.avg_perr;
for fi = 1:size(SAP)
    subplot(6,10,fi)
    sap = reshape(squeeze(SAP(fi,:,:)),length(rp{3}), length(rp{1}),[]);
    contourf(X,Y,sap(:,:,noise_ind)./sap(1,1,noise_ind),20,'LineWidth',1)
    xlabel('lc')
    ylabel('ld')
%     title(tit{fi})
    colormap("jet")
    colorbar;
end

end