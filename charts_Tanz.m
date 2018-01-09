CM = colormap(copper(length(prev_HR(1,:))+3));

figure(1)
hold on; grid on;
for i = 1:length(prev_HR(1,:))
    h(i) = plot(prev_HR(:,i), 'color', CM(i+2,:), 'linewidth',2);
end
h2 = plot(prev_HR(:,2)', 'color', 'b', 'linewidth',2);
set(gca, 'Ylim',[0,90],'YTick',0:10:90,'YTickLabel',0:10:90, 'Xlim',[3,23],'XTick',3:5:23,'XTickLabel',2010:5:2030);
ylabel('HCV prevalence among PWID (%)');
title({['Impact of harm reduction on projected HCV prevalence among PWID' 10 '\rm\fontsize{10}Scaled up over three years and maintained']});
legend([h(1);h2;h(3:end)],'0% coverage','6% coverage (estimated current)','10% coverage','20%','30%','40%','50% coverage');
hold off;

figure(2)
hold on; grid on;
for i = 1:length(inc_HR(1,:))
    h(i) = plot(inc_HR(:,i), 'color', CM(i+2,:), 'linewidth',2);
end
h2 = plot(inc_HR(:,2)', 'color', 'b', 'linewidth',2);
set(gca, 'Ylim',[0,10000],'YTick',0:2500:10000,'YTickLabel',0:2500:10000, 'Xlim',[3,23],'XTick',3:5:23,'XTickLabel',2010:5:2030);
ylabel('HCV incidence');
title({['Impact of harm reduction on projected HCV incidence' 10 '\rm\fontsize{10}Scaled up over three years and maintained']});
legend([h(1);h2;h(3:end)],'0% coverage','6% coverage (estimated current)','10% coverage','15%','20%','25%','30%','35%','40%','45%','50% coverage');
hold off;

