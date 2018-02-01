CM = colormap(copper(length(prev_HR(1,:))+3));

figure(1)
hold on; grid on;
for i = 1:length(prev_HR(1,:))
    h(i) = plot(prev_HR(:,i), 'color', CM(i+2,:), 'linewidth',2);
end
h2 = plot(prev_HR(:,2)', 'color', 'b', 'linewidth',2);
d1 = scatter(8,22,100,'filled','d','k');
set(gca, 'Ylim',[0,90],'YTick',0:10:90,'YTickLabel',0:10:90, 'Xlim',[3,23],'XTick',3:5:23,'XTickLabel',2010:5:2030);
ylabel('HCV prevalence among PWID (%)');
title({['Impact of harm reduction on projected HCV prevalence among PWID' 10 '\rm\fontsize{10}Scaled up over three years and maintained']});
legend([d1;h(1);h2;h(3:end)],'Data','0% coverage','6% coverage (estimated current)','10% coverage','20%','30%','40%','50% coverage');
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
legend([h(1);h2;h(3:end)],'0% coverage','6% coverage (estimated current)','10% coverage','20%','30%','40%','50% coverage');
hold off;


err_ub = abs((0*[inc_HR(:,:)]));% ...
    %- inc_HR(:,:)));
err_lb = abs((-0*[inc_HR(:,:)]));% ...
    %+ inc_HR(:,:)));
figure(3)
% subplot(1,2,1)
set(gca,'LineStyleOrder',{'o','*','h','^','d','s','p'},'ColorOrder', [0 0 0],'NextPlot','replacechildren','XTick',1:5:100,'XTickLabel',2010:5:2050,'FontSize',12);
hold on; grid on;
% for i = 1:length(inc_HR(1,:))
% if i ~=2 
%     h(i)=errorbar(1:length(inc_HR(:,1))',inc_HR(:,i)',err_lb(:,i)',err_ub(:,i)', 'MarkerSize',8); 
% end
% end
% h(2)=errorbar(1:length(inc_HR(:,1))',inc_HR(:,2)',err_lb(:,2)',err_ub(:,2)', 'MarkerSize',8); 
for i = 1:length(inc_HR(1,:))
if i ~=2 
    h(i)=plot([-10*inc_HR(1:10,i);inc_HR(11:end,i)], 'Markersize', 8); 
end
end
h(2)=plot(inc_HR(:,2)','Markersize',8); 
set(h(2),'MarkerEdgecolor','b');
%h = plot(inc_year_per', 'MarkerSize',8); hold on;
%h2 = plot(repmat((1-target_inc)*inc_year_per(1,6),length(inc_year_per)),'--k');
set(gca, 'Ylim',[0,6000],'YTick',0:2000:10000,'YTickLabel',0:2000:10000, 'Xlim',[3,23],'XTick',3:5:23,'XTickLabel',2010:5:2030);
ylabel('HCV incidence');
title({['Impact of harm reduction on projected HCV incidence' 10 '\rm\fontsize{10}Scaled up over three years and maintained']});
legend([h(1);h(2);h(3:end)],'0% coverage','6% coverage (estimated current)','10% coverage','20%','30%','40%','50% coverage','location','northwest');
hold off;




%% LIVER DISEASE AND CARE CASCADE PROJECTION AND BURN IN / DISEASE DISTRIBUTION CALIBRATION RESULTS
scen1=[sum(sum(sum(ycomb_noage(:,:,:,1:5),2),3),4),sum(sum(ycomb_noage(:,:,:,6),2),3),sum(sum(sum(ycomb_noage(:,:,:,7:11),2),3),4),reshape(sum(sum(ycomb_noage(:,:,:,12:18),2),3),length(ycomb_noage(:,1,1,1)),7),sum(sum(sum(ycomb_noage(:,:,:,19:20),2),3),4)];
prev_PWID_scen1=100*sum(sum(ycomb_noage(:,1,:,6:20),3),4)./sum(sum(ycomb_noage(:,1,:,1:20),3),4);
scen1_cas=[sum(sum(ycomb_noage(:,:,1,6:20),2),4),sum(sum(ycomb_noage(:,:,2,6:20),2),4),sum(sum(ycomb_noage(:,:,3,6:20),2),4),sum(sum(ycomb_noage(:,:,4,6:20),2),4),sum(sum(ycomb_noage(:,:,5,6:20),2),4),sum(sum(ycomb_noage(:,:,6,6:20),2),4),sum(sum(ycomb_noage(:,:,7,6:20),2),4),sum(sum(ycomb_noage(:,:,8,6:20),2),4),sum(sum(ycomb_noage(:,:,9,6:20),2),4),sum(sum(ycomb_noage(:,:,10,6:20),2),4)];
 
dat = disease0(2:end)/sum(disease0(2:end));
dat = [dat(1:3),sum(dat(4:5))]';
year=0;
casc = zeros(1,4);
casc(1) = sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,6:13))))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,6:20))));
casc(2) = sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,14)))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,6:20))));
casc(3) = sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,15)))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,6:20))));
casc(4) = sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,16:20))))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,6:20))));

mod = casc';

figure(20)
subplot(1,2,1)
grid on; hold on;
a1_=area(TT2_treat,sum(scen1(:,4:11),2),'Facecolor',CM(8,:));
a2_=area(TT2_treat,sum(scen1(:,5:11),2),'Facecolor',CM(7,:));
a3_=area(TT2_treat,sum(scen1(:,6:11),2),'Facecolor',CM(6,:));
a4_=area(TT2_treat,sum(scen1(:,7:11),2),'Facecolor',CM(5,:));
a5_=area(TT2_treat,sum(scen1(:,8:11),2),'Facecolor',CM(4,:));
a6_=area(TT2_treat,sum(scen1(:,9:11),2),'Facecolor',CM(3,:));
a7_=area(TT2_treat,sum(scen1(:,10:11),2),'Facecolor',CM(2,:));
a8_=area(TT2_treat,sum(scen1(:,11:11),2),'Facecolor',CM(1,:));
d1 = scatter(65,60000,100,'filled','d','k');
xlabel('\fontsize{14}Year'); ylabel({'\fontsize{14} Chronic HCV infections (thousand)'},'Color','k'); 
%[extra,a9__,a9_]=plotyy(TT2_treat,prev_PWID_scen1,TT2_treat,prev_PWID_scen1); 
%set(a9_,'Color','k','LineStyle','--','Linewidth',2); set(a9__,'visible','off');
%set(extra(2), 'XTick',[],'xlim',[50,80],'ylim',[0,65],'YTick',0:10:65,'ycolor','k','FontSize',12);
set(gca, 'xlim',[50,80],'ylim',[0,70000],'XTick',0:10:100,'XTickLabel',1950:10:2050,'FontSize',12,'YTick',0:10000:100000,'YTickLabel',0:10:100,'ycolor','k');
%set(get(extra(2),'Ylabel'),'String','Prevalence among PWID (%)','fontsize',14);
leg=legend([d1,a1_,a2_,a3_,a4_,a5_,a6_,a7_,a8_],{'Data','F0','F1','F2','F3','F4','DC','HCC','LT'},'location','Northwest');
%legpatch=findobj(leg,'type','patch'); hatchfill(legpatch,'single',45,1,'blue');
title({'\fontsize{12}Projected people living with HCV,','by liver disease stage'}); hold off;
subplot(1,2,2)
b = bar(1:2,fliplr([mod,dat]'),'stacked'); hold on;
title({'Model calibration against', '2015 fibrosis data'}, 'fontsize', 12);
ylabel({'Fibrosis stage of people with chronic HCV'});
legend([b(4),b(3),b(2),b(1)],'F0-F1','F2','F3', '>=F4', 'Location', 'Northeast')
set(gca,'XTickLabel',{'Model', 'Data'},'fontsize',12,...
    'YTick',0:.10:1,'YTickLabel',0:10:100)
ylim([0,1]);
b(1).FaceColor=CM(3,:); b(2).FaceColor=CM(5,:); b(3).FaceColor=CM(6,:); b(4).FaceColor=CM(7,:);
axes('Position',[0 0 1 1],'Visible','off');
text(0.5,0.98,{'\fontsize{16}Baseline projections for HCV infections and liver disease in Dar es Salaam'},'HorizontalAlignment','Center')
hold off;