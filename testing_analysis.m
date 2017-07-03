clear all
global filename scenario sens target_late Tin Run start_year num_pops num_cascade num_age num_intervention num_engagement num_region  ... %meta parameters
    mu mu_PWID mu_former exit_IDU r_relapse delta alpha alpha_old alpha_DAA  omega age_cohort P y0 t0...
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2... % disease progression parameters
    Q_sens q_svr q_svr_PWID q_treat  c_daa discount treat  ... % QALY and cost parameters
    p_complete total_PWID PWID0 age_mix prev0 ...%PWID parameters
    incarceration_rate prison_length ... %Prison parameters
    cascade0 cascade0_PWID disease0 cases0 ost0 nsp0 HCC0 diagnoses0 infect_factor infect_base progression_base imported... % data and calibration
    ost_duration  nsp_duration treat_projected target_inc target_death cascade_scale_time care RNAtesting... %intervention
    infect progression imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9 ost_enrollment nsp_enrollment % calibtation parameters


load('C:\Users\Nick\Desktop\Matlab Sims\Testing\sample')

infect_base=infect;
progression_base = progression;

%Run model in prior to 2016
[TT1,y1]=DE_track_age(Tin,y0,t0,treat);
y1_end=reshape(y1(end,:,:,:,:,:,:,:), num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region,33);
death2015 = (sum(sum(sum(sum(sum(sum(y1(end,:,:,:,:,:,:,22))))))) - sum(sum(sum(sum(sum(sum(y1(find(TT1>=65,1),:,:,:,:,:,:,22)))))))) / ...
    (TT1(end)-TT1(find(TT1>=65,1)));
inc2015 = sum(sum(sum(sum(sum(sum(sum(y1(find(TT1>=65,1)+1:end,:,:,:,:,:,:,27))))))))/(66 - TT1(find(TT1>=65,1)+1));
y1_end(:,:,:,:,:,:,21:(27+6))=0;
y1_end(:,6,:,:,:,:,1:20) = y1_end(:,6,:,:,:,:,1:20) + y1_end(:,8,:,:,:,:,1:20) + y1_end(:,10,:,:,:,:,1:20);
y1_end(:,8,:,:,:,:,1:20) = 0; y1_end(:,10,:,:,:,:,1:20) = 0; % moving failed to be re-eligible with DAAs

%% Followup rates for RNA testing after Ab testing frequencies
scenario = 'current';
progression = progression_base;
progression(1,5,2,1) = 1/(30/365); progression(2,5,2,1) = 1/(30/365); progression(3,5,2,1) = 1/(30/365);
progression(1,6,2,1) = 1/(30/365); progression(2,6,2,1) = 1/(30/365); progression(3,6,2,1) = 1/(30/365);
progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
target_late = -1; 
alpha = alpha_DAA;
range_test = [-1, 0.5, 1, 2];
range_followup = [-1,0:0.2:1];
treat = treat_projected;
for i = 1:length(range_test)
    for j = 1:length(range_followup)
        if range_test(i) > 0
            temp_rate = progression_base(1,1,2,1);
            progression(1,1,2,1) = range_test(i);% + 0.1*temp_rate;
            progression(2,1,2,1) = range_test(i);% + 0.1*temp_rate;
            %progression(3,1,2,1) = 0.9*range_test(i) + 0.1*temp_rate;
        end
        if range_followup(j) > 0
            temp_rate1 = (progression_base(1,2,2,1)  - 0.46*(1/0.25))/(1-0.46); % 46% followed up in 3 months (0.25 years) currently;
            temp_rate2 = (progression_base(2,2,2,1)  - 0.46*(1/0.25))/(1-0.46); % 46% followed up in 3 months (0.25 years) currently
            progression(1,2,2,1) = max(0,range_followup(j)*(1/0.25) + (1-range_followup(j))*temp_rate1);
            progression(2,2,2,1) = max(0,range_followup(j)*(1/0.25) + (1-range_followup(j))*temp_rate2);
        else
            progression(1,2,2,1) = progression_base(1,2,2,1);
            progression(2,2,2,1) = progression_base(2,2,2,1);
        end
        
        [TT2_treat3,y2_treat3]=DE_track_age(Run,y1_end,TT1,treat);
        [ycomb3_noage, summary(3,:,s), tr3, tr3_] = gather_outputs(y1,y2_treat3,TT2_treat3);
        
        for year_elim = 1:Run-2
            incred_test(i,year_elim,j) = (inc2015 - sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin+year_elim,1)+1:find(TT2_treat3>=Tin+year_elim+1,1),:,:,27)))))/...
                (TT2_treat3(find(TT2_treat3>=Tin+year_elim+1,1))-TT2_treat3(find(TT2_treat3>=Tin+year_elim,1))) / inc2015;
        end
        if incred_test(i,end) >= 0.8
            year_elim_test(i,j) = 1950+Tin+find(incred_test(i,:,j)>=0.8,1);
        else
            year_elim_test(i,j) = -1;
        end
    end
end


%%  RNA testing
scenarios = {'RNA','RNA_all'};
for j = 1:length(scenarios)
    scenario = scenarios(j);
    progression = progression_base;
    progression(1,5,2,1) = 1/(30/365); progression(2,5,2,1) = 1/(30/365); progression(3,5,2,1) = 1/(30/365);
    progression(1,6,2,1) = 1/(30/365); progression(2,6,2,1) = 1/(30/365); progression(3,6,2,1) = 1/(30/365);
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    target_late = -1;
    treat = treat_projected;
    
    RNAtesting_range = [-1,1/2,1,2];
    for i = 1:length(RNAtesting_range)
        if RNAtesting_range(i) > 0
            RNAtesting = RNAtesting_range(i);
        else
            temp_rate1 = progression_base(1,2,2,1);
            RNAtesting = temp_rate1;
        end
        [TT2_treat4,y2_treat4]=DE_track_age(Run,y1_end,TT1,treat);
        [ycomb4_noage, summary(4,:,s), tr4, tr4_] = gather_outputs(y1,y2_treat4,TT2_treat4);
        for year_elim = 1:Run-2
            incred_RNA(i,year_elim,j) = (inc2015 - sum(sum(sum(ycomb4_noage(find(TT2_treat4>=Tin+year_elim,1)+1:find(TT2_treat4>=Tin+year_elim+1,1),:,:,27)))))/...
                (TT2_treat4(find(TT2_treat4>=Tin+year_elim+1,1))-TT2_treat4(find(TT2_treat4>=Tin+year_elim,1))) / inc2015;
        end
        if incred_RNA(i,end,j) >= 0.8
            year_elim_RNA(i,j) = 1950+Tin+find(incred_RNA(i,:,j)>=0.8,1);
        else
            year_elim_RNA(i) = -1;
        end
    end
end

%% Plots
a0 = 100*reshape(incred_test(:,end,2:end),length(range_test),length(range_followup)-1);
for i = 1:length(a0(:,1))
    for j = 2:length(a0(1,:))
        a(i,j) = max(0, min(a0(i,j)-a0(i,1), a0(i,j) - a0(i,j-1)));
    end
end
a1 = 100*reshape(incred_RNA(:,end,:),length(RNAtesting_range),length(scenarios));
for i = 1:length(a1(:,1))
    for j = 1:length(a1(1,:))
        a2(i,j) = a1(i,j) - sum(a1(i,1:j-1));
    end
end

CM1 = colormap(copper(9));
CM2 =colormap(summer(9));
CM = [CM1;CM2];
figure(1)
b = bar(a,'stacked'); grid on;
b(1).Parent.Parent.Colormap = flipud(CM(3:8,:));
set(b,'edgecolor','none');
ylim([0,100]); ylabel('Estimated incidence reduction in 2030 (%)');
legend('Current average RNA follow-up rates', '20% RNA tested within 3 months', '40% RNA tested within 3 months', '60% RNA tested within 3 months', '80% RNA tested within 3 months', '100% RNA tested within 3 months','location','southeast');
set(gca,'XTicklabel',{'Current \newlinetesting rates','Two-yearly \newlineAb testing, \newlineall PWID',...
    'Annual \newlineAb testing, \newlineall PWID',...
    'Six-monthly \newlineAb testing, \newlineall PWID'});
title(['\fontsize{14}2030 incidence reduction following treatment availability' char(10) '\it{}\fontsize{12}Changes to antibody testing and follow-up rates among PWID'],'HorizontalAlignment','center');

figure(2)
b2 = bar(a2,'stacked'); grid on;
b2(1).Parent.Parent.Colormap = flipud(CM(5:end,:)); set(b2,'edgecolor','none');
ylim([0,100]); ylabel('Estimated incidence reduction in 2030 (%)');
legend('Coverage: PWID on OST', 'Coverage: all PWID', '40% within 3 months', '60% within 3 months', '80% within 3 months', '100% within 3 months','location','southeast');
set(gca,'XTicklabel',{'Current Ab\newlinetesting rates \newlinereplaced with \newlineRNA tests','Two-yearly \newlineRNA testing', 'Annual \newlineRNA testing', 'Six-monthly \newlineRNA testing'});
title(['\fontsize{14}2030 incidence reduction following treatment availability' char(10) '\it{}\fontsize{12}Antibody testing replaced with RNA testing among PWID'],'HorizontalAlignment','center');

figure(3)
subplot(1,2,1)
b = bar(a,'stacked'); grid on;
b(1).Parent.Parent.Colormap = flipud(CM1(3:8,:));
set(b,'edgecolor','none');
ylim([0,100]); ylabel('Estimated incidence reduction in 2030 (%)');
legend('Current average RNA follow-up rates', '20% RNA tested within 3 months', '40% RNA tested within 3 months', '60% RNA tested within 3 months', '80% RNA tested within 3 months', '100% RNA tested within 3 months','location','southeast');
set(gca,'XTicklabel',{'Current \newlinetesting rates','Two-yearly \newlineAb testing, \newlineall PWID',...
    'Annual \newlineAb testing, \newlineall PWID',...
    'Six-monthly \newlineAb testing, \newlineall PWID'});
title(['\it{}\fontsize{10}Changes to antibody testing and follow-up rates among PWID'],'HorizontalAlignment','center');
hold off;
subplot(1,2,2)
b2 = bar(a2,'stacked'); grid on;
b2(1).Parent.Parent.Colormap = flipud(CM1(3:8,:)); set(b2,'edgecolor','none');
ylim([0,100]); ylabel('Estimated incidence reduction in 2030 (%)');
legend('Coverage: PWID on OST', 'Coverage: all PWID', '40% within 3 months', '60% within 3 months', '80% within 3 months', '100% within 3 months','location','southeast');
set(gca,'XTicklabel',{'Current Ab\newlinetesting rates \newlinereplaced with \newlineRNA tests','Two-yearly \newlineRNA testing', 'Annual \newlineRNA testing', 'Six-monthly \newlineRNA testing'});
title(['\it{}\fontsize{10}Antibody testing replaced with RNA testing among PWID'],'HorizontalAlignment','center');
axes('Position',[0 0 1 1],'Visible','off');
text(0.5,0.98,{'\fontsize{14}2030 incidence reduction following treatment availability'},'HorizontalAlignment','Center')



%% Calibrate to 30% prevalence
% scenario = 'empty';
% loaddata
% prev0 = [2015, 0.3;...
%     2005, 0.4];
% alpha=alpha_old; % calibrate in pre-DAA era
% target_late=1; % target treatments to people with late-liver disease
% 
% [output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases] = calibrate_optim_par(15, 8);
% save('C:\Users\Nick\Desktop\Matlab Sims\Testing\30prevcalibration')
