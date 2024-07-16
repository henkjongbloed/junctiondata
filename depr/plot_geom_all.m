function plot_geom_all(U)
figure;
hold on
for j = 1:numel(U)
    U{j}.xs.plot()
    U{j}.B.plot()
%     plot3(S{j}.Xp', S{j}.Yp', S{j}.Z' , '.')

end
grid on
xlabel('x')
ylabel('y')
hold off

end