function plot_geom_U(U)

subplot(1, 2, 1)
hold on
for b = 1:numel(U)
    U{b}.xs.plot()
    U{b}.B.plot()
    %plot(cell2mat(S{b}.x'), cell2mat(S{b}.y'), '.')    
end
grid on
hold off
subplot(1, 2, 2)
hold on
for b = 1:numel(U)
    plot(U{1, b}.eta)
    BN{b} = U{b}.BN;
end
legend(BN(:))
hold off

end