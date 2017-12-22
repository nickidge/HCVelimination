total_PWID_model=sum(sum(ycomb2_noage(:,1,:,1:20),3),4);
total_former_model=sum(sum(ycomb2_noage(:,2,:,1:20),3),4);
total_import_model=sum(sum(ycomb2_noage(:,3,:,1:20),3),4);
total_PWID_model_inf=sum(sum(ycomb2_noage(:,1,:,6:20),3),4);
total_former_model_inf=sum(sum(ycomb2_noage(:,2,:,6:20),3),4);
total_import_model_inf=sum(sum(ycomb2_noage(:,3,:,6:20),3),4);
total_inf=sum(sum(sum(ycomb2_noage(:,1:3,:,6:20),2),3),4);

CM = colormap(copper(9));
CM2 = colormap(summer(9));

figure(1) % Model population sizes over time
set(gca, 'ColorOrder', [0,0,1;0,0.5,0;1,0,0], 'LineStyleOrder',{'-','--'},...
    'NextPlot', 'replacechildren','XTick',1:10:100,'XTickLabel',1950:10:2050); grid on;
plot(TT2_treat,[total_PWID_model,total_former_model,total_import_model,total_PWID_model_inf,total_former_model_inf,total_import_model_inf]); hold on;
plot(TT2_treat, total_inf,'-k');
legend('Total PWID', 'Total former', 'Total non-PWID','Infected PWID', 'Infected former', 'Infected non-PWID', 'Total infections','location','Northwest');
%ylim([0,250000]);
ylabel('Number in model (1000)');
title('Population sizes in model');

figure(4) %Which sub-populations are treatments actually going to?
p=plot(treatment_alloc'); hold on;
plot(sum(treatment_alloc,1),'r');
set(p(1),'LineStyle','--','Marker','s','Markersize',6,'Color','k','LineWidth',2);
set(p(2),'LineStyle','-','Marker','d','Markersize',6,'Color','k','LineWidth',2);
set(p(3),'LineStyle','--','Marker','s','Markersize',6,'Color',[160,160,160]/255,'LineWidth',2);
set(p(4),'LineStyle','-','Marker','d','Markersize',6,'Color',[160,160,160]/255,'LineWidth',2);
set(gca,'Fontsize',12,'ColorOrder',[0,0,0],'XTick',0:5:15,'XTicklabel',2015:5:2030)
legend('PWID early liver disease','PWID advanced liver disease','Former PWID early liver disease','Former PWID advanced liver disease','location','West')
xlim([1,15]); grid on; hold off;
%ylim([0,6000]); 
title('\fontsize{14}Prioritising to current PWID')
ylabel('Annual number of treatments')
xlabel('Year')

%% CALIBRATION RESULTS FOR CARE CASCADE
dat = cascade0(1,2:end);
casc = zeros(3,4);
casc(1,1) = sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1,2:10,6:20)))./sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1,1:10,6:20)));
casc(1,2) = sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1,3:10,6:20)))./sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1,1:10,6:20)));
casc(1,3) = sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1,4:10,6:20)))./sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1,1:10,6:20)));
casc(1,4) = sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1,6:10,6:20)))./sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1,1:10,6:20)));
casc(2,1) = sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),2:3,2:10,6:20))))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),2:3,1:10,6:20))));
casc(2,2) = sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),2:3,3:10,6:20))))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),2:3,1:10,6:20))));
casc(2,3) = sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),2:3,4:10,6:20))))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),2:3,1:10,6:20))));
casc(2,4) = sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),2:3,6:10,6:20))))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),2:3,1:10,6:20))));
casc(3,1) = sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1:3,2:10,6:20))))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1:3,1:10,6:20))));
casc(3,2) = sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1:3,3:10,6:20))))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1:3,1:10,6:20))));
casc(3,3) = sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1:3,4:10,6:20))))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1:3,1:10,6:20))));
casc(3,4) = sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1:3,6:10,6:20))))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),1:3,1:10,6:20))));


mod = casc';
figure(6)
b = bar(1:4,[mod(1:4,:),dat(1:4)']);
title('Model calibration results for 2015 cascade of care', 'fontsize', 12);
xlabel('Step in cascade');
ylabel({'Percentage of people with chronic HCV', 'who have reached cascade stage'});
legend('Model: PWID','Model: non-IDU','Model: overall', 'Data', 'Location', 'Northeast')
set(gca,'XTickLabel',{'Antibody+', 'RNA+', 'Genotyped', 'Liver test'},...
    'YTick',0:.10:1,'YTickLabel',0:10:100)
ylim([0,1]);
b(1).FaceColor=CM2(6,:); b(2).FaceColor=CM2(4,:); b(3).FaceColor=CM2(3,:); b(4).FaceColor=CM2(1,:);
b(1).EdgeColor=CM2(6,:); b(2).EdgeColor=CM2(4,:); b(3).EdgeColor=CM2(3,:); b(4).EdgeColor=CM2(1,:);
b(1).LineWidth = 1; b(2).LineWidth = 1; b(3).LineWidth = 1; b(3).LineWidth = 1;
text(0.75,mod(1,1),strcat(num2str(100*mod(1,1),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(0.91,mod(1,2),strcat(num2str(100*mod(1,2),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(1.1,mod(1,3),strcat(num2str(100*mod(1,3),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(1.3,dat(1),strcat(num2str(100*dat(1),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(1.75,mod(2,1),strcat(num2str(100*mod(2,1),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(1.91,mod(2,2),strcat(num2str(100*mod(2,2),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(2.1,mod(2,3),strcat(num2str(100*mod(2,3),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(2.3,dat(2),strcat(num2str(100*dat(2),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(2.75,mod(3,1),strcat(num2str(100*mod(3,1),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(2.91,mod(3,2),strcat(num2str(100*mod(3,2),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(3.1,mod(3,3),strcat(num2str(100*mod(3,3),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(3.3,dat(3),strcat(num2str(100*dat(3),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(3.75,mod(4,1),strcat(num2str(100*mod(4,1),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(3.91,mod(4,2),strcat(num2str(100*mod(4,2),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(4.1,mod(4,3),strcat(num2str(100*mod(4,3),2),'%'),'horiz','center','vert','bottom','Fontsize',8)
text(4.3,dat(4),strcat(num2str(100*dat(4),2),'%'),'horiz','center','vert','bottom','Fontsize',8)

%% DISEASE DISTRIBUTION CALIBRATION RESULTS
dat = disease0(2:end)/sum(disease0(2:end));
dat = [dat(1:3),sum(dat(4:5))]';
year=0;
casc = zeros(1,4);
casc(1) = sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,6:13))))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,6:20))));
casc(2) = sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,14)))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,6:20))));
casc(3) = sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,15)))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,6:20))));
casc(4) = sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,16:20))))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + year,1),1:3,1:10,6:20))));

mod = casc';
figure(8)
b = bar(1:2,fliplr([mod,dat]'),'stacked');
title('Model calibration results for 2015 fibrosis levels', 'fontsize', 12);
ylabel({'Fibrosis stage of people with chronic HCV'});
legend([b(4),b(3),b(2),b(1)],'F0-F1','F2','F3', '>=F4', 'Location', 'Northeast')
set(gca,'XTickLabel',{'Model', 'Data'},'fontsize',12,...
    'YTick',0:.10:1,'YTickLabel',0:10:100)
ylim([0,1]);
b(1).FaceColor=CM(3,:); b(2).FaceColor=CM(5,:); b(3).FaceColor=CM(6,:); b(4).FaceColor=CM(7,:);



%% WHO targets
inc_year_per = 100+100*(inc_year - repmat(inc_year(:,1),1,length(inc_year(1,:))))./repmat(inc_year(:,1),1,length(inc_year(1,:)));
death_year_per = 100+100*(death_year - repmat(death_year(:,10),1,length(death_year(1,:))))./repmat(death_year(:,10),1,length(death_year(1,:)));

inc_year_per = [inc_year(1,:); inc_year(5,:); inc_year([2,4],:)] * 100 ./ repmat(inc_year([1,5,2,4],6),1,length(inc_year(1,:)));
death_year_per = [death_year(1,:); death_year(5:6,:); death_year(2:3,:)];

figure(5)
% subplot(1,2,1)
set(gca,'LineStyleOrder',{'o','*','h','^','-'},'ColorOrder', [0 0 0],'NextPlot','replacechildren','XTick',1:5:100,'XTickLabel',2010:5:2050,'FontSize',12);
h = plot(inc_year_per', 'MarkerSize',8); hold on;
h2 = plot(repmat((1-target_inc)*inc_year_per(1,6),length(inc_year_per)),'--k');
l = legend([h(1);h(3);h(2);h(4);h2(1)],{'Baseline',...
    ['Scenario 2: projected uptake'],...
    ['Scenario 3: min. PWID treatment requirements'],...    
    ['Scenario 4: projected uptake + annual testing'],...
    ['WHO target']});%,...
set(l,'fontsize',10);
title({['Annual HCV incidence (among PWID)' 10 '\fontsize{10}For various treatment scenarios']})
ylabel('Relative incidence compared to 2015 (%)'); xlabel('Year');
xlim([1,21]);grid on;
% subplot(1,2,2)
% set(gca,'LineStyleOrder',{'o','*','h','^','+','-'},'ColorOrder',[0,0,0],'NextPlot','replacechildren','XTick',1:5:100,'XTickLabel',2010:5:2050,'FontSize',12);
% h = plot(death_year_per', 'Markersize', 8); hold on;
% h2 = plot(repmat((1-target_death)*death_year_per(1,6),length(inc_year_per)),'--k');
% title({['Annual liver-related mortality']})
% ylabel('Annual liver-related deaths'); xlabel('Year');
% xlim([1,21]);  grid on;
%axes('Position',[0 0 1 1],'Visible','off');
%text(0.5,0.98,{'\fontsize{16}HCV incidence (among PWID) and liver-related deaths for various cascade programs'},'HorizontalAlignment','Center')




%% LIVER DISEASE AND CARE CASCADE PROJECTION AND BURN IN
scen1=[sum(sum(sum(ycomb_noage(:,:,:,1:5),2),3),4),sum(sum(ycomb_noage(:,:,:,6),2),3),sum(sum(sum(ycomb_noage(:,:,:,7:11),2),3),4),reshape(sum(sum(ycomb_noage(:,:,:,12:18),2),3),length(ycomb_noage(:,1,1,1)),7),sum(sum(sum(ycomb_noage(:,:,:,19:20),2),3),4)];
prev_PWID_scen1=100*sum(sum(ycomb_noage(:,1,:,6:20),3),4)./sum(sum(ycomb_noage(:,1,:,1:20),3),4);
scen1_cas=[sum(sum(ycomb_noage(:,:,1,6:20),2),4),sum(sum(ycomb_noage(:,:,2,6:20),2),4),sum(sum(ycomb_noage(:,:,3,6:20),2),4),sum(sum(ycomb_noage(:,:,4,6:20),2),4),sum(sum(ycomb_noage(:,:,5,6:20),2),4),sum(sum(ycomb_noage(:,:,6,6:20),2),4),sum(sum(ycomb_noage(:,:,7,6:20),2),4),sum(sum(ycomb_noage(:,:,8,6:20),2),4),sum(sum(ycomb_noage(:,:,9,6:20),2),4),sum(sum(ycomb_noage(:,:,10,6:20),2),4)];

figure(10)
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
xlabel('\fontsize{14}Year'); ylabel({'\fontsize{14} Chronic HCV infections (thousand)'},'Color','k'); 
[extra,a9__,a9_]=plotyy(TT2_treat,prev_PWID_scen1,TT2_treat,prev_PWID_scen1); 
set(a9_,'Color','k','LineStyle','--','Linewidth',2); set(a9__,'visible','off');
set(extra(2), 'XTick',[],'xlim',[0,80],'ylim',[0,65],'YTick',0:10:65,'ycolor','k','FontSize',12);
set(extra(1), 'xlim',[0,80],'ylim',[0,1200000],'XTick',0:10:100,'XTickLabel',1950:10:2050,'FontSize',12,'YTick',0:100000:1200000,'YTickLabel',0:100:1200,'ycolor','k');
set(get(extra(2),'Ylabel'),'String','Prevalence among PWID (%)','fontsize',14);
leg=legend([a1_,a2_,a3_,a4_,a5_,a6_,a7_,a8_,a9_],{'F0','F1','F2','F3','F4','DC','HCC','LT',['Prevalence among' 10 'current PWID']},'location','Northwest');
%legpatch=findobj(leg,'type','patch'); hatchfill(legpatch,'single',45,1,'blue');
title({'\fontsize{12}Liver disease stages'}); hold off;
subplot(1,2,2)
grid on; hold on;
a1_=area(TT2_treat,sum(scen1_cas(:,1:10),2),'Facecolor',CM2(7,:));
a2_=area(TT2_treat,sum(scen1_cas(:,2:10),2),'Facecolor',CM2(6,:));
a3_=area(TT2_treat,sum(scen1_cas(:,3:10),2),'Facecolor',CM2(5,:));
a4_=area(TT2_treat,sum(scen1_cas(:,4:10),2),'Facecolor',CM2(4,:));
a5_=area(TT2_treat,sum(scen1_cas(:,6:10),2),'Facecolor',CM2(3,:));
a6_=area(TT2_treat,sum(scen1_cas(:,7:10),2),'Facecolor',CM2(2,:));
a7_=area(TT2_treat,sum(scen1_cas(:,[8,10]),2),'Facecolor',CM2(1,:));
xlabel('\fontsize{14}Year'); %ylabel({'\fontsize{14}IDU-acquired infections in liver', 'disease stage (thousand)'},'Color','k'); 
[extra,a9__,a9_]=plotyy(TT2_treat,prev_PWID_scen1,TT2_treat,prev_PWID_scen1); 
set(a9_,'Color','k','LineStyle','--','Linewidth',2); set(a9__,'visible','off');
set(extra(2), 'XTick',[],'xlim',[0,80],'ylim',[0,65],'YTick',0:10:60,'ycolor','k','FontSize',12);
set(extra(1), 'xlim',[0,80],'ylim',[0,1200000],'XTick',0:10:100,'XTickLabel',1950:10:2050,'FontSize',12,'YTick',0:100000:1200000,'YTickLabel',0:100:1200,'ycolor','k');
set(get(extra(2),'Ylabel'),'String','Prevalence among PWID (%)','fontsize',14);
%legpatch=findobj(leg,'type','patch'); hatchfill(legpatch,'single',45,1,'blue');
title({'\fontsize{12}Scenario 1: basic community-based care'}); hold off;
legend([a1_,a2_,a3_,a4_,a5_,a6_,a7_,a9_],{'Undiagnosed','HCV antibody+','HCV RNA+','Genotyped','Liver disease tested','Commenced treatment','Failed treatment',['Prevalence among' 10 'current PWID']},'location','Southeast')
title('\fontsize{12}Cascade stages')
axes('Position',[0 0 1 1],'Visible','off');
text(0.5,0.98,{'\fontsize{16}Australian HCV infections','\fontsize{14}Projected outcomes to scale-up to 90% diagnosed by 2030'},'HorizontalAlignment','Center')
hold off;


%% Disease burden
base=[sum(sum(sum(ycomb_noage(:,:,:,1:5),2),3),4),sum(sum(ycomb_noage(:,:,:,6),2),3),sum(sum(sum(ycomb_noage(:,:,:,7:11),2),3),4),reshape(sum(sum(ycomb_noage(:,:,:,12:18),2),3),length(ycomb_noage(:,1,1,1)),7),sum(sum(sum(ycomb_noage(:,:,:,19:20),2),3),4)];
scen1=[sum(sum(sum(ycomb2_noage(:,:,:,1:5),2),3),4),sum(sum(ycomb2_noage(:,:,:,6),2),3),sum(sum(sum(ycomb2_noage(:,:,:,7:11),2),3),4),reshape(sum(sum(ycomb2_noage(:,:,:,12:18),2),3),length(ycomb2_noage(:,1,1,1)),7),sum(sum(sum(ycomb2_noage(:,:,:,19:20),2),3),4)];
scen2=[sum(sum(sum(ycomb3_noage(:,:,:,1:5),2),3),4),sum(sum(ycomb3_noage(:,:,:,6),2),3),sum(sum(sum(ycomb3_noage(:,:,:,7:11),2),3),4),reshape(sum(sum(ycomb3_noage(:,:,:,12:18),2),3),length(ycomb3_noage(:,1,1,1)),7),sum(sum(sum(ycomb3_noage(:,:,:,19:20),2),3),4)];
scen3=[sum(sum(sum(ycomb4_noage(:,:,:,1:5),2),3),4),sum(sum(ycomb4_noage(:,:,:,6),2),3),sum(sum(sum(ycomb4_noage(:,:,:,7:11),2),3),4),reshape(sum(sum(ycomb4_noage(:,:,:,12:18),2),3),length(ycomb4_noage(:,1,1,1)),7),sum(sum(sum(ycomb4_noage(:,:,:,19:20),2),3),4)];
scen4=[sum(sum(sum(ycomb2_noage(:,:,:,1:5),2),3),4),sum(sum(ycomb2_noage(:,:,:,6),2),3),sum(sum(sum(ycomb2_noage(:,:,:,7:11),2),3),4),reshape(sum(sum(ycomb2_noage(:,:,:,12:18),2),3),length(ycomb2_noage(:,1,1,1)),7),sum(sum(sum(ycomb2_noage(:,:,:,19:20),2),3),4)];

prev_PWID_base=100*sum(sum(ycomb_noage(:,1,:,6:20),3),4)./sum(sum(ycomb_noage(:,1,:,1:20),3),4);
prev_PWID_scen1=100*sum(sum(ycomb2_noage(:,1,:,6:20),3),4)./sum(sum(ycomb2_noage(:,1,:,1:20),3),4);
prev_PWID_scen2=100*sum(sum(ycomb3_noage(:,1,:,6:20),3),4)./sum(sum(ycomb3_noage(:,1,:,1:20),3),4);
prev_PWID_scen3=100*sum(sum(ycomb4_noage(:,1,:,6:20),3),4)./sum(sum(ycomb4_noage(:,1,:,1:20),3),4);
prev_PWID_scen5=100*sum(sum(ycomb5_noage(:,1,:,6:20),3),4)./sum(sum(ycomb5_noage(:,1,:,1:20),3),4);
prev_PWID_scen6=100*sum(sum(ycomb6_noage(:,1,:,6:20),3),4)./sum(sum(ycomb6_noage(:,1,:,1:20),3),4);

TT2_treatopt_min=TT2_treat;

figure(2)
subplot(1,2,1)
grid on; hold on;
%a0=area(TT2,sum(base(:,1:11),2),'Facecolor',[255,255,255]/255); 
a1=area(TT2,sum(base(:,4:11),2),'Facecolor',CM(8,:));
a2=area(TT2,sum(base(:,5:11),2),'Facecolor',CM(7,:));
a3=area(TT2,sum(base(:,6:11),2),'Facecolor',CM(6,:));
a4=area(TT2,sum(base(:,7:11),2),'Facecolor',CM(5,:));
a5=area(TT2,sum(base(:,8:11),2),'Facecolor',CM(4,:));
a6=area(TT2,sum(base(:,9:11),2),'Facecolor',CM(3,:));
a7=area(TT2,sum(base(:,10:11),2),'Facecolor',CM(2,:));
a8=area(TT2,sum(base(:,11:11),2),'Facecolor',CM(1,:));
title('\fontsize{12}No treatment available')
xlabel('\fontsize{14}Year'); ylabel({'\fontsize{14}HCV infections in liver disease stage'},'Color','k'); 
[extra,a9_,a9]=plotyy(TT2,prev_PWID_base,TT2,prev_PWID_base); 
set(a9,'Color','k','LineStyle','--','Linewidth',2); set(a9_,'visible','off');
set(extra(2), 'XTick',[],'xlim',[50,80],'ylim',[0,65],'YTick',0:10:65,'ycolor','k','FontSize',12);
set(extra(1), 'xlim',[50,80],'ylim',[0,250000],'XTick',50:10:100,'XTickLabel',2000:10:2050,'FontSize',12,'YTick',0:50000:250000,'YTickLabel',0:50:250,'ycolor','k');
%set(get(extra(2),'Ylabel'),'String','Prevalence among PWID (%)','fontsize',14);
legend([a1,a2,a3,a4,a5,a6,a7,a8,a9],{'F0','F1','F2','F3','F4','DC','HCC','LT',['Prevalence among' 10 'current PWID']},'location','Northwest')
subplot(1,2,2)
grid on; hold on;
a1_=area(TT2_treat,sum(scen1(:,4:11),2),'Facecolor',CM(8,:));
a2_=area(TT2_treat,sum(scen1(:,5:11),2),'Facecolor',CM(7,:));
a3_=area(TT2_treat,sum(scen1(:,6:11),2),'Facecolor',CM(6,:));
a4_=area(TT2_treat,sum(scen1(:,7:11),2),'Facecolor',CM(5,:));
a5_=area(TT2_treat,sum(scen1(:,8:11),2),'Facecolor',CM(4,:));
a6_=area(TT2_treat,sum(scen1(:,9:11),2),'Facecolor',CM(3,:));
a7_=area(TT2_treat,sum(scen1(:,10:11),2),'Facecolor',CM(2,:));
a8_=area(TT2_treat,sum(scen1(:,11:11),2),'Facecolor',CM(1,:));
xlabel('\fontsize{14}Year'); %ylabel({'\fontsize{14}IDU-acquired infections in liver', 'disease stage (thousand)'},'Color','k'); 
[extra,a9__,a9_]=plotyy(TT2_treat,prev_PWID_scen1,TT2_treat,prev_PWID_scen1); 
set(a9_,'Color','k','LineStyle','--','Linewidth',2); set(a9__,'visible','off');
set(extra(2), 'XTick',[],'xlim',[50,80],'ylim',[0,65],'YTick',0:10:65,'ycolor','k','FontSize',12);
set(extra(1), 'xlim',[50,80],'ylim',[0,250000],'XTick',50:10:100,'XTickLabel',2000:10:2050,'FontSize',12,'YTick',0:50000:250000,'YTickLabel',0:50:250,'ycolor','k');
set(get(extra(2),'Ylabel'),'String','Prevalence among PWID (%)','fontsize',14);
%leg=legend([a1_,a2_,a3_,a4_,a5_,a6_,a7_,a8_,a9_],{'F0','F1','F2','F3','F4','DC','HCC','LT',['Prevalence among' 10 'current PWID']},'location','Northwest');
%legpatch=findobj(leg,'type','patch'); hatchfill(legpatch,'single',45,1,'blue');
title({'\fontsize{12}Treat PWID to hit WHO incidence target'})
axes('Position',[0 0 1 1],'Visible','off');
text(0.5,0.98,{'\fontsize{16}HCV related liver disease in Iceland','\fontsize{14}Projected outcomes 2015-2030'},'HorizontalAlignment','Center')



%% CARE CASCADE UNDER DIFFERENT SCENARIOS
base=[sum(sum(ycomb_noage(:,:,1,6:20),2),4),sum(sum(ycomb_noage(:,:,2,6:20),2),4),sum(sum(ycomb_noage(:,:,3,6:20),2),4),sum(sum(ycomb_noage(:,:,4,6:20),2),4),sum(sum(ycomb_noage(:,:,5,6:20),2),4),sum(sum(ycomb_noage(:,:,6,6:20),2),4),sum(sum(ycomb_noage(:,:,7,6:20),2),4),sum(sum(ycomb_noage(:,:,8,6:20),2),4),sum(sum(ycomb_noage(:,:,9,6:20),2),4),sum(sum(ycomb_noage(:,:,10,6:20),2),4)];
scen1=[sum(sum(ycomb2_noage(:,:,1,6:20),2),4),sum(sum(ycomb2_noage(:,:,2,6:20),2),4),sum(sum(ycomb2_noage(:,:,3,6:20),2),4),sum(sum(ycomb2_noage(:,:,4,6:20),2),4),sum(sum(ycomb2_noage(:,:,5,6:20),2),4),sum(sum(ycomb2_noage(:,:,6,6:20),2),4),sum(sum(ycomb2_noage(:,:,7,6:20),2),4),sum(sum(ycomb2_noage(:,:,8,6:20),2),4),sum(sum(ycomb2_noage(:,:,9,6:20),2),4),sum(sum(ycomb2_noage(:,:,10,6:20),2),4)];
scen2=[sum(sum(ycomb3_noage(:,:,1,6:20),2),4),sum(sum(ycomb3_noage(:,:,2,6:20),2),4),sum(sum(ycomb3_noage(:,:,3,6:20),2),4),sum(sum(ycomb3_noage(:,:,4,6:20),2),4),sum(sum(ycomb3_noage(:,:,5,6:20),2),4),sum(sum(ycomb3_noage(:,:,6,6:20),2),4),sum(sum(ycomb3_noage(:,:,7,6:20),2),4),sum(sum(ycomb3_noage(:,:,8,6:20),2),4),sum(sum(ycomb3_noage(:,:,9,6:20),2),4),sum(sum(ycomb3_noage(:,:,10,6:20),2),4)];
scen3=[sum(sum(ycomb4_noage(:,:,1,6:20),2),4),sum(sum(ycomb4_noage(:,:,2,6:20),2),4),sum(sum(ycomb4_noage(:,:,3,6:20),2),4),sum(sum(ycomb4_noage(:,:,4,6:20),2),4),sum(sum(ycomb4_noage(:,:,5,6:20),2),4),sum(sum(ycomb4_noage(:,:,6,6:20),2),4),sum(sum(ycomb4_noage(:,:,7,6:20),2),4),sum(sum(ycomb4_noage(:,:,8,6:20),2),4),sum(sum(ycomb4_noage(:,:,9,6:20),2),4),sum(sum(ycomb4_noage(:,:,10,6:20),2),4)];
scen5=[sum(sum(ycomb5_noage(:,:,1,6:20),2),4),sum(sum(ycomb5_noage(:,:,2,6:20),2),4),sum(sum(ycomb5_noage(:,:,3,6:20),2),4),sum(sum(ycomb5_noage(:,:,4,6:20),2),4),sum(sum(ycomb5_noage(:,:,5,6:20),2),4),sum(sum(ycomb5_noage(:,:,6,6:20),2),4),sum(sum(ycomb5_noage(:,:,7,6:20),2),4),sum(sum(ycomb5_noage(:,:,8,6:20),2),4),sum(sum(ycomb5_noage(:,:,9,6:20),2),4),sum(sum(ycomb5_noage(:,:,10,6:20),2),4)];
scen6=[sum(sum(ycomb6_noage(:,:,1,6:20),2),4),sum(sum(ycomb6_noage(:,:,2,6:20),2),4),sum(sum(ycomb6_noage(:,:,3,6:20),2),4),sum(sum(ycomb6_noage(:,:,4,6:20),2),4),sum(sum(ycomb6_noage(:,:,5,6:20),2),4),sum(sum(ycomb6_noage(:,:,6,6:20),2),4),sum(sum(ycomb6_noage(:,:,7,6:20),2),4),sum(sum(ycomb6_noage(:,:,8,6:20),2),4),sum(sum(ycomb6_noage(:,:,9,6:20),2),4),sum(sum(ycomb6_noage(:,:,10,6:20),2),4)];

  
figure(3)
subplot(2,2,1)
grid on; grid minor; hold on;
%a0=area(TT2,sum(base(:,1:11),2),'Facecolor',[255,255,255]/255); 
a1=area(TT2,sum(base(:,1:10),2),'Facecolor',[224,224,224]/255);
a2=area(TT2,sum(base(:,2:10),2),'Facecolor',[192,192,192]/255);
a3=area(TT2,sum(base(:,3:10),2),'Facecolor',[160,160,160]/255);
a4=area(TT2,sum(base(:,4:10),2),'Facecolor',[128,128,128]/255);
a5=area(TT2,sum(base(:,6:10),2),'Facecolor',[96,96,96]/255);
a6=area(TT2,sum(base(:,7:10),2),'Facecolor',[64,64,64]/255);
a7=area(TT2,sum(base(:,[8,10]),2),'Facecolor',[32,32,32]/255);
%a8=area(TT2,sum(base(:,8:8),2),'Facecolor',[0,0,0]/255);
title('\fontsize{12}Baseline')
%xlabel('\fontsize{14}Year'); 
ylabel('\fontsize{14}Infections in cascade stage (thousand)','Color','k'); 
[extra,a9_,a9]=plotyy(TT2,prev_PWID_base,TT2,prev_PWID_base); 
set(a9,'Color','k','LineStyle','--','Linewidth',2); set(a9_,'visible','off');
set(extra(2), 'XTick',[],'xlim',[60,80],'ylim',[0,60],'YTick',0:10:60,'ycolor','k','FontSize',12);
set(extra(1), 'xlim',[60,80],'ylim',[0,250000],'XTick',0:5:100,'XTickLabel',1950:5:2050,'FontSize',12,'YTick',0:25000:250000,'YTickLabel',0:25:250,'ycolor','k');
%set(get(extra(2),'Ylabel'),'String','Prevalence among PWID (%)','fontsize',14);
%legend([a1,a2,a3,a4,a5,a6,a7,a9],{'Undiagnosed','HCV antibody+','HCV RNA+','Genotyped','Liver disease tested','Commenced treatment','Failed treatment',['Prevalence among' 10 'current PWID']},'location','Northwest')
hold off;
subplot(2,2,2)
grid on; grid minor; hold on;
a1_=area(TT2_treat5,sum(scen5(:,1:10),2),'Facecolor',[224,224,224]/255);
a2_=area(TT2_treat5,sum(scen5(:,2:10),2),'Facecolor',[192,192,192]/255);
a3_=area(TT2_treat5,sum(scen5(:,3:10),2),'Facecolor',[160,160,160]/255);
a4_=area(TT2_treat5,sum(scen5(:,4:10),2),'Facecolor',[128,128,128]/255);
a5_=area(TT2_treat5,sum(scen5(:,6:10),2),'Facecolor',[96,96,96]/255);
a6_=area(TT2_treat5,sum(scen5(:,7:10),2),'Facecolor',[64,64,64]/255);
a7_=area(TT2_treat5,sum(scen5(:,[8,10]),2),'Facecolor',[32,32,32]/255);
xlabel('\fontsize{14}Year'); %ylabel({'\fontsize{14}IDU-acquired infections in liver', 'disease stage (thousand)'},'Color','k'); 
[extra,a9__,a9_]=plotyy(TT2_treat5,prev_PWID_scen5,TT2_treat5,prev_PWID_scen5); 
set(a9_,'Color','k','LineStyle','--','Linewidth',2); set(a9__,'visible','off');
set(extra(2), 'XTick',[],'xlim',[60,80],'ylim',[0,60],'YTick',0:10:60,'ycolor','k','FontSize',12);
set(extra(1), 'xlim',[60,80],'ylim',[0,250000],'XTick',0:5:100,'XTickLabel',1950:5:2050,'FontSize',12,'YTick',0:25000:250000,'YTickLabel',0:25:250,'ycolor','k');
set(get(extra(2),'Ylabel'),'String','Prevalence among PWID (%)','fontsize',14);
%legpatch=findobj(leg,'type','patch'); hatchfill(legpatch,'single',45,1,'blue');
title({'\fontsize{12}Scenario 2: scaled up primary care + APRI'}); hold off;
subplot(2,2,3)
grid on; grid minor; hold on;
a1_F2=area(TT2_treat6,sum(scen6(:,1:10),2),'Facecolor',[224,224,224]/255);
a2_F2=area(TT2_treat6,sum(scen6(:,2:10),2),'Facecolor',[192,192,192]/255);
a3_F2=area(TT2_treat6,sum(scen6(:,3:10),2),'Facecolor',[160,160,160]/255);
a4_F2=area(TT2_treat6,sum(scen6(:,4:10),2),'Facecolor',[128,128,128]/255);
a5_F2=area(TT2_treat6,sum(scen6(:,6:10),2),'Facecolor',[96,96,96]/255);
a6_F2=area(TT2_treat6,sum(scen6(:,7:10),2),'Facecolor',[64,64,64]/255);
a7_F2=area(TT2_treat6,sum(scen6(:,[8,10]),2),'Facecolor',[32,32,32]/255);
%xlabel('\fontsize{14}Year'); 
%ylabel('\fontsize{14}Infections in cascade stage (thousand)','Color','k'); 
[extra,a9_F2_,a9_F2]=plotyy(TT2_treat6,prev_PWID_scen6,TT2_treat6,prev_PWID_scen6); 
set(a9_F2,'Color','k','LineStyle','--','Linewidth',2); set(a9_F2_,'visible','off');
set(extra(2), 'XTick',[],'xlim',[60,80],'ylim',[0,60],'YTick',0:10:60,'ycolor','k','FontSize',12);
set(extra(1), 'xlim',[60,80],'ylim',[0,250000],'XTick',0:5:100,'XTickLabel',1950:5:2050,'FontSize',12,'YTick',0:25000:250000,'YTickLabel',0:25:250,'ycolor','k');
%set(get(extra(2),'Ylabel'),'String','Prevalence among PWID (%)','fontsize',14);
title({'\fontsize{12}Scenario 3: scaled up primary care + APRI', '+ annual testing of PWID on OST'}); hold off;
subplot(2,2,4)
grid on; grid minor; hold on;
a1_3=area(TT2_treat4,sum(scen3(:,1:10),2),'Facecolor',[224,224,224]/255);
a2_3=area(TT2_treat4,sum(scen3(:,2:10),2),'Facecolor',[192,192,192]/255);
a3_3=area(TT2_treat4,sum(scen3(:,3:10),2),'Facecolor',[160,160,160]/255);
a4_3=area(TT2_treat4,sum(scen3(:,4:10),2),'Facecolor',[128,128,128]/255);
a5_3=area(TT2_treat4,sum(scen3(:,6:10),2),'Facecolor',[96,96,96]/255);
a6_3=area(TT2_treat4,sum(scen3(:,7:10),2),'Facecolor',[64,64,64]/255);
a7_3=area(TT2_treat4,sum(scen3(:,[8,10]),2),'Facecolor',[32,32,32]/255);
xlabel('\fontsize{14}Year'); %ylabel({'\fontsize{14}IDU-acquired infections in liver', 'disease stage (thousand)'},'Color','k'); 
[extra,a9_3_,a9_3]=plotyy(TT2_treat4,prev_PWID_scen3,TT2_treat4,prev_PWID_scen3); 
set(a9_3,'Color','k','LineStyle','--','Linewidth',2); set(a9_3_,'visible','off');
legend([a1,a2,a3,a4,a5,a6,a7,a9],{'Undiagnosed','HCV antibody+','HCV RNA+','Genotyped','Liver disease tested','Commenced treatment','Failed treatment',['Prevalence among' 10 'current PWID']},'location','Southeast')

set(extra(2), 'XTick',[],'xlim',[60,80],'ylim',[0,60],'YTick',0:10:60,'ycolor','k','FontSize',12);
set(extra(1), 'xlim',[60,80],'ylim',[0,250000],'XTick',0:5:100,'XTickLabel',1950:5:2050,'FontSize',12,'YTick',0:25000:250000,'YTickLabel',0:25:250,'ycolor','k');
set(get(extra(2),'Ylabel'),'String','Prevalence among PWID (%)','fontsize',14);
title({'\fontsize{12}Scenario 5: all health system interventions'})
axes('Position',[0 0 1 1],'Visible','off');
text(0.5,0.98,{'\fontsize{16}Australian HCV infections in various cascade stages','\fontsize{14}Projected outcomes 2016-2030 under different scenarios'},'HorizontalAlignment','Center')
hold off;



%% UNUSED
% %% HCC and diagnosis calibration
% HCC_year_per = [HCC_year(1:2,:); HCC_year(5:6,:); HCC_year(3:4,:)];
% diagnosed_year_per = [diagnosed_year(1:2,:); diagnosed_year(5:6,:); diagnosed_year(3:4,:)];
% HCC_data = [0,0,0,0,0,flipud(HCC0(1:9,2))']';
% diagnosed_data = flipud(diagnoses0(1:10,2));
% 
% figure(11)
% subplot(1,2,1)
% h2 =  plot(HCC_data,'*r'); hold on;
% h = plot(HCC_year_per(2,:)','r'); 
% h_base = plot(HCC_year_per(1,:)','--r'); 
% title({['Annual HCC incidence']})
% ylabel('New cases of HCV-related HCC in year'); xlabel('Year');
% xlim([6,26]); ylim([0,10]); grid on;
% set(gca,'XTick',1:5:100,'XTickLabel',2000:5:2050,'FontSize',12); hold off;
% legend([h2,h],'Data','Model');
% subplot(1,2,2)
% h2 =  plot(diagnosed_data,'*b'); hold on;
% h = plot(diagnosed_year_per(2,:),'b'); 
% title({['Annual HCV RNA+ diagnoses']})
% ylabel('Number of diagnoses in year'); xlabel('Year');
% xlim([6,26]); ylim([0,100]); grid on;
% set(gca,'XTick',1:5:100,'XTickLabel',2000:5:2050,'FontSize',12); hold off;
% legend([h2,h],'Data','Model');



% 
% figure(9)
% grid on; hold on;
% %a0=area(TT2,sum(base(:,1:11),2),'Facecolor',[255,255,255]/255); 
% a1=area(TT2,sum(base(:,4:11),2),'Facecolor',CM(8,:));
% a2=area(TT2,sum(base(:,5:11),2),'Facecolor',CM(7,:));
% a3=area(TT2,sum(base(:,6:11),2),'Facecolor',CM(6,:));
% a4=area(TT2,sum(base(:,7:11),2),'Facecolor',CM(5,:));
% a5=area(TT2,sum(base(:,8:11),2),'Facecolor',CM(4,:));
% a6=area(TT2,sum(base(:,9:11),2),'Facecolor',CM(3,:));
% a7=area(TT2,sum(base(:,10:11),2),'Facecolor',CM(2,:));
% a8=area(TT2,sum(base(:,11:11),2),'Facecolor',CM(1,:));
% %title('\fontsize{12}No treatment available')
% xlabel('\fontsize{14}Year'); ylabel({'\fontsize{14}HCV infections in liver disease stage'},'Color','k'); 
% [extra,a9_,a9]=plotyy(TT2,prev_PWID_base,TT2,prev_PWID_base); 
% set(a9,'Color','k','LineStyle','--','Linewidth',2); set(a9_,'visible','off');
% set(extra(2), 'XTick',[],'xlim',[40,66],'ylim',[0,65],'YTick',0:10:60,'ycolor','k','FontSize',12);
% set(extra(1), 'xlim',[40,66],'ylim',[0,1750],'XTick',0:5:100,'XTickLabel',1950:5:2050,'FontSize',12,'YTick',0:250:2500,'YTickLabel',0:250:2500,'ycolor','k');
% set(get(extra(2),'Ylabel'),'String','Prevalence among PWID (%)','fontsize',14);
% legend([a1,a2,a3,a4,a5,a6,a7,a8,a9],{'F0','F1','F2','F3','F4','DC','HCC','LT',['Prevalence among' 10 'current PWID']},'location','Northwest')
% axes('Position',[0 0 1 1],'Visible','off');
% text(0.5,0.98,{'\fontsize{16}Modelled HCV related liver disease in Iceland','\fontsize{14}Model burn-in with no DAA treatments available'},'HorizontalAlignment','Center')
% hold off;
% 






% %% Projected cascade
% year=15;
% for t=1:year+1
%     g(1,t,1) = sum(sum(ycomb_noage(find(TT2>=Tin + t-1,1),:,1,6:20)))./sum(sum(sum(ycomb_noage(find(TT2>=Tin,1),:,1:10,6:20))));
%     g(2,t,1) = sum(sum(ycomb_noage(find(TT2>=Tin + t-1,1),:,2,6:20)))./sum(sum(sum(ycomb_noage(find(TT2>=Tin,1),:,1:10,6:20))));
%     g(3,t,1) = sum(sum(ycomb_noage(find(TT2>=Tin + t-1,1),:,3,6:20)))./sum(sum(sum(ycomb_noage(find(TT2>=Tin,1),:,1:10,6:20))));
%     g(4,t,1) = sum(sum(ycomb_noage(find(TT2>=Tin + t-1,1),:,4,6:20)))./sum(sum(sum(ycomb_noage(find(TT2>=Tin,1),:,1:10,6:20))));
%     g(5,t,1) = sum(sum(ycomb_noage(find(TT2>=Tin + t-1,1),:,6,6:20)))./sum(sum(sum(ycomb_noage(find(TT2>=Tin,1),:,1:10,6:20))));
%     g(6,t,1) = sum(sum(sum(ycomb_noage(find(TT2>=Tin + t-1,1),:,[7,8,9],6:20))))./sum(sum(sum(ycomb_noage(find(TT2>=Tin,1),:,1:10,6:20))));
%     g(8,t,1) = sum(death_year(1,1:t))./sum(sum(sum(ycomb_noage(find(TT2>=Tin,1),:,1:10,6:20))));
%     g(7,t,1) = 1-sum(g(1:8,t,1));
%     
%     g(1,t,2) = sum(sum(ycomb2_noage(find(TT2_treat>=Tin + t-1,1),:,1,6:20)))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),:,1:10,6:20))));
%     g(2,t,2) = sum(sum(ycomb2_noage(find(TT2_treat>=Tin + t-1,1),:,2,6:20)))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),:,1:10,6:20))));
%     g(3,t,2) = sum(sum(ycomb2_noage(find(TT2_treat>=Tin + t-1,1),:,3,6:20)))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),:,1:10,6:20))));
%     g(4,t,2) = sum(sum(ycomb2_noage(find(TT2_treat>=Tin + t-1,1),:,4,6:20)))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),:,1:10,6:20))));
%     g(5,t,2) = sum(sum(ycomb2_noage(find(TT2_treat>=Tin + t-1,1),:,6,6:20)))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),:,1:10,6:20))));
%     g(6,t,2) = sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin + t-1,1),:,[7,8,9],6:20))))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),:,1:10,6:20))));
%     g(8,t,2) = sum(death_year(2,1:t))./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),:,1:10,6:20))));
%     g(7,t,2) = 1-sum(g(1:8,t,2));% 1000*total_treat(4,2)./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),:,1:10,6:20))));
%     
%     g(1,t,3) = sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + t-1,1),:,1,6:20)))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin,1),:,1:10,6:20))));
%     g(2,t,3) = sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + t-1,1),:,2,6:20)))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin,1),:,1:10,6:20))));
%     g(3,t,3) = sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + t-1,1),:,3,6:20)))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin,1),:,1:10,6:20))));
%     g(4,t,3) = sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + t-1,1),:,4,6:20)))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin,1),:,1:10,6:20))));
%     g(5,t,3) = sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + t-1,1),:,6,6:20)))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin,1),:,1:10,6:20))));
%     g(6,t,3) = sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin + t-1,1),:,[7,8,9],6:20))))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin,1),:,1:10,6:20))));
%     g(8,t,3) = sum(death_year(3,1:t))./sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin,1),:,1:10,6:20))));
%     g(7,t,3) = 1-sum(g(1:8,t,3));% 1000*total_treat(4,2)./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),:,1:10,6:20))));
%     
%     g(1,t,4) = sum(sum(ycomb4_noage(find(TT2_treat4>=Tin + t-1,1),:,1,6:20)))./sum(sum(sum(ycomb4_noage(find(TT2_treat4>=Tin,1),:,1:10,6:20))));
%     g(2,t,4) = sum(sum(ycomb4_noage(find(TT2_treat4>=Tin + t-1,1),:,2,6:20)))./sum(sum(sum(ycomb4_noage(find(TT2_treat4>=Tin,1),:,1:10,6:20))));
%     g(3,t,4) = sum(sum(ycomb4_noage(find(TT2_treat4>=Tin + t-1,1),:,3,6:20)))./sum(sum(sum(ycomb4_noage(find(TT2_treat4>=Tin,1),:,1:10,6:20))));
%     g(4,t,4) = sum(sum(ycomb4_noage(find(TT2_treat4>=Tin + t-1,1),:,4,6:20)))./sum(sum(sum(ycomb4_noage(find(TT2_treat4>=Tin,1),:,1:10,6:20))));
%     g(5,t,4) = sum(sum(ycomb4_noage(find(TT2_treat4>=Tin + t-1,1),:,6,6:20)))./sum(sum(sum(ycomb4_noage(find(TT2_treat4>=Tin,1),:,1:10,6:20))));
%     g(6,t,4) = sum(sum(sum(ycomb4_noage(find(TT2_treat4>=Tin + t-1,1),:,[7,8,9],6:20))))./sum(sum(sum(ycomb4_noage(find(TT2_treat4>=Tin,1),:,1:10,6:20))));
%     g(8,t,4) = sum(death_year(4,1:t))./sum(sum(sum(ycomb4_noage(find(TT2_treat4>=Tin,1),:,1:10,6:20))));
%     g(7,t,4) = 1-sum(g(1:8,t,4));% 1000*total_treat(4,2)./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),:,1:10,6:20))));
%     
%     g(1,t,5) = sum(sum(ycomb5_noage(find(TT2_treat5>=Tin + t-1,1),:,1,6:20)))./sum(sum(sum(ycomb5_noage(find(TT2_treat5>=Tin,1),:,1:10,6:20))));
%     g(2,t,5) = sum(sum(ycomb5_noage(find(TT2_treat5>=Tin + t-1,1),:,2,6:20)))./sum(sum(sum(ycomb5_noage(find(TT2_treat5>=Tin,1),:,1:10,6:20))));
%     g(3,t,5) = sum(sum(ycomb5_noage(find(TT2_treat5>=Tin + t-1,1),:,3,6:20)))./sum(sum(sum(ycomb5_noage(find(TT2_treat5>=Tin,1),:,1:10,6:20))));
%     g(4,t,5) = sum(sum(ycomb5_noage(find(TT2_treat5>=Tin + t-1,1),:,4,6:20)))./sum(sum(sum(ycomb5_noage(find(TT2_treat5>=Tin,1),:,1:10,6:20))));
%     g(5,t,5) = sum(sum(ycomb5_noage(find(TT2_treat5>=Tin + t-1,1),:,6,6:20)))./sum(sum(sum(ycomb5_noage(find(TT2_treat5>=Tin,1),:,1:10,6:20))));
%     g(6,t,5) = sum(sum(sum(ycomb5_noage(find(TT2_treat5>=Tin + t-1,1),:,[7,8,9],6:20))))./sum(sum(sum(ycomb5_noage(find(TT2_treat5>=Tin,1),:,1:10,6:20))));
%     g(8,t,5) = sum(death_year(5,1:t))./sum(sum(sum(ycomb5_noage(find(TT2_treat5>=Tin,1),:,1:10,6:20))));
%     g(7,t,5) = 1-sum(g(1:8,t,5));% 1000*total_treat(4,2)./sum(sum(sum(ycomb2_noage(find(TT2_treat>=Tin,1),:,1:10,6:20))));
% end
% 
% 
% figure(7)
% b = bar(fliplr(g(:,:,5)'),'stacked');
% title({'\fontsize{14}Projected cascade of care','\fontsize{12}Scaled up primary care + APRI + annual testing of PWID on OST'},'FontWeight','Normal');
% xlabel('Year');
% ylabel({'Percentage of people who had chronic HCV in 2015'});
% l = {'Cumulative deaths','Cured','Treated','Liver tested and treatment ready','Genotyped / needing a liver assessment','RNA+ diagnosed needing a genotype test','Antibody+ needing an RNA test','Undiagnosed'};
% legend(fliplr(b), fliplr(l),'Location', 'East','fontsize',12)
% set(gca,'XTick',1:5:16,'XTickLabel',2015:5:2030, ... %{},...
%     'YTick',0:.10:1,'YTickLabel',0:10:100,'fontsize',12)
% ylim([0,1]); xlim([0.5,16.5]);
% b(8).FaceColor = [224,224,224]/255;
% b(7).FaceColor = [192,192,192]/255;
% b(6).FaceColor = [160,160,160]/255;
% b(5).FaceColor = [128,128,128]/255;
% b(4).FaceColor = [96,96,96]/255;
% b(3).FaceColor = [64,64,64]/255;
% b(2).FaceColor = [0,0,0]/255;
% b(1).FaceColor = [255,255,255]/255;
% b(1).FaceColor = [0,0,0];
% b(2).FaceColor = [0,0.7,0];
% b(3).FaceColor = [0,0.8,0.8];
% b(4).FaceColor = [0,0.5,0.8];
% b(5).FaceColor = [0.5,0.5,0.8];
% b(6).FaceColor = [1,0.9,0];
% b(7).FaceColor = [1,0.6,0];
% b(8).FaceColor = [1,0,0];
% 
% figure(7)
% subplot(2,2,1)
% b = bar(fliplr(g(:,:,1)'),'stacked');
% title('Baseline', 'fontsize', 12);
% xlabel('Year');
% ylabel({'Percentage of people who had chronic HCV in 2015'});
% l = {'Died','Cured','Treated','Liver tested waiting for treatment','Genotyped needing a liver test','RNA+ diagnosed needing a genotype test','Antibody+ needing an RNA test','Undiagnosed'};
% legend(fliplr(b), fliplr(l),'Location', 'East')
% set(gca,'XTick',1:5:16,'XTickLabel',2015:5:2030, ... %{},...
%     'YTick',0:.10:1,'YTickLabel',0:10:100)
% ylim([0,1]); xlim([0.5,16.5]);
% b(1).FaceColor = [0,0,0];
% b(2).FaceColor = [0,0.7,0];
% b(3).FaceColor = [0,0.8,0.8];
% b(4).FaceColor = [0,0.5,0.8];
% b(5).FaceColor = [0.5,0.5,0.8];
% b(6).FaceColor = [1,0.9,0];
% b(7).FaceColor = [1,0.6,0];
% b(8).FaceColor = [1,0,0];
% subplot(2,2,2)
% b = bar(fliplr(g(:,:,5)'),'stacked');
% title('Scaled up primary care + APRI', 'fontsize', 12);
% xlabel('Year');
% ylabel({'Percentage of people who had chronic HCV in 2015'});
% l = {'Died','Cured','Treated','Liver tested waiting for treatment','Genotyped needing a liver test','RNA+ diagnosed needing a genotype test','Antibody+ needing an RNA test','Undiagnosed'};
% legend(fliplr(b), fliplr(l),'Location', 'East')
% set(gca,'XTick',1:5:16,'XTickLabel',2015:5:2030, ... %{},...
%     'YTick',0:.10:1,'YTickLabel',0:10:100)
% ylim([0,1]); xlim([0.5,16.5]);
% b(1).FaceColor = [0,0,0];
% b(2).FaceColor = [0,0.7,0];
% b(3).FaceColor = [0,0.8,0.8];
% b(4).FaceColor = [0,0.5,0.8];
% b(5).FaceColor = [0.5,0.5,0.8];
% b(6).FaceColor = [1,0.9,0];
% b(7).FaceColor = [1,0.6,0];
% b(8).FaceColor = [1,0,0];
% subplot(2,2,3)
% b = bar(fliplr(g(:,:,3)'),'stacked');
% title('Scaled up primary care + APRI + rapid RNA', 'fontsize', 12);
% xlabel('Year');
% ylabel({'Percentage of people who had chronic HCV in 2015'});
% l = {'Died','Cured','Treated','Liver tested waiting for treatment','Genotyped needing a liver test','RNA+ diagnosed needing a genotype test','Antibody+ needing an RNA test','Undiagnosed'};
% legend(fliplr(b), fliplr(l),'Location', 'East')
% set(gca,'XTick',1:5:16,'XTickLabel',2015:5:2030, ... %{},...
%     'YTick',0:.10:1,'YTickLabel',0:10:100)
% ylim([0,1]); xlim([0.5,16.5]);
% b(1).FaceColor = [0,0,0];
% b(2).FaceColor = [0,0.7,0];
% b(3).FaceColor = [0,0.8,0.8];
% b(4).FaceColor = [0,0.5,0.8];
% b(5).FaceColor = [0.5,0.5,0.8];
% b(6).FaceColor = [1,0.9,0];
% b(7).FaceColor = [1,0.6,0];
% b(8).FaceColor = [1,0,0];
% subplot(2,2,4)
% b = bar(fliplr(g(:,:,4)'),'stacked');
% title('Scaled up primary care + APRI + rapid RNA + annual testing PWID', 'fontsize', 12);
% xlabel('Year');
% ylabel({'Percentage of people who had chronic HCV in 2015'});
% l = {'Died','Cured','Treated','Liver tested waiting for treatment','Genotyped needing a liver test','RNA+ diagnosed needing a genotype test','Antibody+ needing an RNA test','Undiagnosed'};
% legend(fliplr(b), fliplr(l),'Location', 'East')
% set(gca,'XTick',1:5:16,'XTickLabel',2015:5:2030, ... %{},...
%     'YTick',0:.10:1,'YTickLabel',0:10:100)
% ylim([0,1]); xlim([0.5,16.5]);
% b(1).FaceColor = [0,0,0];
% b(2).FaceColor = [0,0.7,0];
% b(3).FaceColor = [0,0.8,0.8];
% b(4).FaceColor = [0,0.5,0.8];
% b(5).FaceColor = [0.5,0.5,0.8];
% b(6).FaceColor = [1,0.9,0];
% b(7).FaceColor = [1,0.6,0];
% b(8).FaceColor = [1,0,0];


