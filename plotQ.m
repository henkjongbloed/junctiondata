function plotQ(msh)

figure;
subplot(1,2,1)
hold on
plot(msh{1, 1}.cs.Q  )
plot(-msh{2, 1}.cs.Q  )
plot(msh{3, 1}.cs.Q   )

hold off
grid on
xlabel('Time [h]')
ylabel('Q [m^3/s')
title('Water flux through transects')
legend('NM','OM','RWW')

subplot(1,2,2)
hold on
plot(msh{1, 1}.sec.QS )
plot(-msh{2, 1}.sec.QS  )
plot(msh{3, 1}.sec.QS )
hold off
grid on
xlabel('Time [h]')
ylabel('Q_s [kg/s')
title('Salt flux through transects')
legend('NM','OM','RWW')