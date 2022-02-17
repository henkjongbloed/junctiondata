function plot19(TA)

XTL = {'F_{00}', 'F_{01}^*', 'F_{02}^*', 'F_{03}^*', 'F_{04}^*',...
    'F_{10}', 'F_{11}^*',...
    'F_{20}', 'F_{21}^*',...
    'F_{30}', 'F_{31}^*', 'F_{32}^*', 'F_{33}^*', 'F_{34}^*', 'F_{35}', 'F_{36}^*',...
    'F_{40}', 'F_{41}', 'F_{42}'};

BN = {' NM', ' OM', ' NWW'};
br = length(TA);
figure;
for b = 1:br
    
    subplot(2*br,1,2*b-1)
    bar(TA{b}.FT(:,1))
    title(strcat('Lateral Decomposition us at ', BN{b}))
    xticks(1:length(TA{b}.FT(:,1)))
    xticklabels(XTL)
    subplot(2*br,1,2*b)

    bar(TA{b}.FT(:,2))
    title(strcat('Lateral Decomposition vs at ', BN{b}))
    xticks(1:length(TA{b}.FT(:,2)))
    xticklabels(XTL)
end