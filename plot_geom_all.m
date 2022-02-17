function plot_geom_all(U, S, h)

subplot(1, 2, 1)
hold on
for b = 1:numel(U)
    U{b}.xs.plot()
    U{b}.B.plot()
    plot(cell2mat(S{b}.x'), cell2mat(S{b}.y'), '.')    
end
grid on
xlabel('x')
ylabel('y')
hold off
subplot(1, 2, 2)
hold on
plot(h.date, h.level)
st = {'-', '--', '-.'};
for i=1:numel(U)
    plot([min(datenum(U{1, i}.T.adcp.time)), min(datenum(U{1, i}.T.adcp.time))], [min(h.level), max(h.level)], LineStyle= st{i} )
    plot([max(datenum(U{1, i}.T.adcp.time)), max(datenum(U{1, i}.T.adcp.time))], [min(h.level), max(h.level)] , LineStyle= st{i}  )
end
xlabel('time')
ylabel('waterlevel')
hold off

end