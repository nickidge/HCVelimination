clear all
global filename scenario sens target_late Tin Run start_year num_pops num_cascade num_age num_intervention num_engagement num_region  dt ... %meta parameters
    mu mu_PWID mu_former exit_IDU r_relapse delta alpha alpha_old alpha_DAA  omega age_cohort P y0 t0...
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2... % disease progression parameters
    Q_sens q_svr q_svr_PWID q_treat  c_daa discount treat  ... % QALY and cost parameters
    p_complete total_PWID PWID0 age_mix prev0 harm_reduction r_inc_up followup ...%PWID parameters
    incarceration_rate prison_length ... %Prison parameters
    cascade0 cascade0_PWID disease0 cases0 ost0 nsp0 HCC0 diagnoses0 infect_factor infect_base progression_base imported... % data and calibration
    ost_duration  nsp_duration treat_projected target_inc target_death cascade_scale_time care RNAtesting... %intervention
    infect progression imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9 ost_enrollment nsp_enrollment % calibtation parameters


%load('C:\Users\Nick\Desktop\Matlab Sims\Testing\cal_vary_rel_start')
%load('C:\Users\Nick\Desktop\Matlab Sims\Testing\60prev')

load('calibration_data')
load('60prev')
%load('30prev')
dt = 1/12;
s=1;
followup = 1;
infect_base=infect;
progression_base = progression;

%Run model in prior to 2016
[TT1,y1]=DE_track_age(Tin,y0,t0,treat);
y1_end=reshape(y1(end,:,:,:,:,:,:,:), num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region,33);
death2015 = (sum(sum(sum(sum(sum(sum(y1(end,:,:,:,:,:,:,22))))))) - sum(sum(sum(sum(sum(sum(y1(find(TT1>=65,1),:,:,:,:,:,:,22)))))))) / ...
    (TT1(end)-TT1(find(TT1>=65,1)));
inc2015 = sum(sum(sum(sum(sum(sum(sum(y1(find(TT1>=65,1):end-1,:,:,:,:,:,:,27))))))));
%sum(sum(sum(sum(sum(sum(sum(y1(find(TT1>=65,1)+1:end,:,:,:,:,:,:,27))))))))/(66 - TT1(find(TT1>=65,1)+1));
y1_end(:,:,:,:,:,:,21:(27+6))=0;
y1_end(:,6,:,:,:,:,1:20) = y1_end(:,6,:,:,:,:,1:20) + y1_end(:,8,:,:,:,:,1:20) + y1_end(:,10,:,:,:,:,1:20);
y1_end(:,8,:,:,:,:,1:20) = 0; y1_end(:,10,:,:,:,:,1:20) = 0; % moving failed to be re-eligible with DAAs


%% Followup rates for RNA testing after Ab testing frequencies
scenario = 'OST_test';
progression = progression_base;
progression(1,5,2,1) = 1/(30/365); progression(2,5,2,1) = 1/(30/365); progression(3,5,2,1) = 1/(30/365);
progression(1,6,2,1) = 1/(30/365); progression(2,6,2,1) = 1/(30/365); progression(3,6,2,1) = 1/(30/365);
progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
target_late = 1; 
alpha = alpha_DAA;
range_test = [-1, 0.5, 1, 2, 4]; % divided by 2 to assume that infection happens midway between tests
range_followup = [-1,0.46,0.5:0.1:1];
treat = [20000,000,2000];
for i = 1:length(range_test)
    for j = 1:length(range_followup)
        if range_test(i) > 0
            temp_rate = progression_base(1,1,2,1);
            temp_rate0 = progression_base(2,1,2,1);
            progression(1,1,2,1) = range_test(i)*0.9 + 0.1*temp_rate;
            progression(2,1,2,1) = range_test(i)*0.9 + 0.1*temp_rate0;
            %progression(3,1,2,1) = 0.9*range_test(i) + 0.1*temp_rate;
        else
            progression(1,1,2,1) = progression_base(1,1,2,1);
            progression(2,1,2,1) = progression_base(2,1,2,1);
            followup = 1;
        end
        if range_followup(j) > 0
            temp_rate1 = (progression_base(1,2,2,1)  - 0.46*(1/0.25))/(1-0.46); % 46% followed up in 3 months (0.25 years) currently;
            temp_rate2 = (progression_base(2,2,2,1)  - 0.46*(1/0.25))/(1-0.46); % 46% followed up in 3 months (0.25 years) currently
            progression(2,2,2,1) = 1/0.25; progression(1,2,2,1) = 1 / 0.25;
            followup = range_followup(j);%/0.46;

        else
            progression(1,2,2,1) = progression_base(1,2,2,1);
            progression(2,2,2,1) = progression_base(2,2,2,1);
            followup = 1;
        end
        [TT2_treat3,y2_treat3]=DE_track_age(Run,y1_end,TT1,treat);
        [ycomb3_noage, summary(3,:,s), tr3, tr3_] = gather_outputs(y1,y2_treat3,TT2_treat3);
        
        for year_elim = 1:Run-2
            incred_test(i,year_elim,j) = (inc2015 - sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin+year_elim,1):find(TT2_treat3>=Tin+year_elim+1,1)-1,:,:,27)))))/ inc2015;%...
               % (TT2_treat3(find(TT2_treat3>=Tin+year_elim+1,1))-TT2_treat3(find(TT2_treat3>=Tin+year_elim,1))) / inc2015;
        end
        if incred_test(i,end,j) >= 0.8
            year_elim_test(i,j) = 1950+Tin+find(incred_test(i,:,j)>=0.8,1);
        else
            year_elim_test(i,j) = -1;
        end
    end
end


%%  RNA testing
scenarios = {'OST_test','current'};
for j = 1:length(scenarios)
    scenario = scenarios(j);
    progression = progression_base;
    progression(1,5,2,1) = 1/(30/365); progression(2,5,2,1) = 1/(30/365); progression(3,5,2,1) = 1/(30/365);
    progression(1,6,2,1) = 1/(30/365); progression(2,6,2,1) = 1/(30/365); progression(3,6,2,1) = 1/(30/365);
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    target_late = 1;
    treat = [20000,000,2000];
    followup = 1;
    
    RNAtesting_range = [-1,1/2,1,2,4]; % divide by two to represent being infected midway between tests
    for i = 1:length(RNAtesting_range)
        if RNAtesting_range(i) > 0
            progression(2,2,2,1) = 1/0.25; progression(1,2,2,1) = 1 / 0.25;
            progression(1,1,2,1) = range_test(i)*0.9 + 0.1*temp_rate;
            progression(2,1,2,1) = range_test(i)*0.9 + 0.1*temp_rate0;
            followup = 1;
            
            %temp_rate3 = progression_base(1,2,2,1);
            %RNAtesting = RNAtesting_range(i)*0.9 + 0.1*temp_rate3;
            
        else
            temp_rate3 = progression_base(1,2,2,1);
            RNAtesting = temp_rate3;
        end
        [TT2_treat4,y2_treat4]=DE_track_age(Run,y1_end,TT1,treat);
        [ycomb4_noage, summary(4,:,s), tr4, tr4_] = gather_outputs(y1,y2_treat4,TT2_treat4);
        for year_elim = 1:Run-2
            incred_RNA(i,year_elim,j) = (inc2015 - sum(sum(sum(ycomb4_noage(find(TT2_treat4>=Tin+year_elim,1):find(TT2_treat4>=Tin+year_elim+1,1)-1,:,:,27)))))/ inc2015;%...
               % (TT2_treat4(find(TT2_treat4>=Tin+year_elim+1,1))-TT2_treat4(find(TT2_treat4>=Tin+year_elim,1))) / inc2015;
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
a3 = 0*a0;
a3(:,1) = a0(:,1);
a3(:,2) = a0(:,2) - a0(:,1);
a3(:,3) = a0(:,3) - a0(:,2);
a3(:,4) = max(0, a0(:,4) - a0(:,3));
a3(:,5) = max(0, a0(:,5) - a0(:,4));
a3(:,6) = max(0, a0(:,6) - a0(:,5));
a3(:,7) = max(0, a0(:,7) - a0(:,6));


% AVHEC
% a3(:,1) = 100* [0.45, 0.51, 0.62, 0.69];
% a3(:,2) = 6.2 * ones(1,4);
% a3(:,3) = 4.4 * ones(1,4);
% a3(:,4) = 4.4 * ones(1,4);
% a3(:,5) = 4.4 * ones(1,4);
% a3(:,6) = 4.4 * ones(1,4);
% 
% 
% a2(:,1) = 100 * [0.64, 0.74, 0.85, 0.87];
% a2(:,2) = [5,5,4,3];

    
CM1 = colormap(copper(9));
CM2 =colormap(summer(9));
CM = [CM1;CM2];
figure(1)
b = bar(a3,'stacked'); grid on;
b(1).Parent.Parent.Colormap = flipud(CM(3:8,:));
set(b,'edgecolor','none');
ylim([0,100]); ylabel('Estimated incidence reduction in 2030 (%)');
legend('Current average RNA follow-up rates (46%)', '50% RNA tested within 3 months', '60% ...', '70% ...', '80% ...',...
    '90% ...','100% RNA tested within 3 months', 'location','southeast');
set(gca,'XTicklabel',{'Current \newlinetesting rates','Two-yearly \newlineAb testing, \newlineall PWID',...
    'Annual \newlineAb testing, \newlineall PWID',...
    'Six-monthly \newlineAb testing, \newlineall PWID'});
title(['\fontsize{14}2030 incidence reduction following treatment availability' char(10) '\it{}\fontsize{12}Changes to antibody testing and follow-up rates among PWID'],'HorizontalAlignment','center');

figure(2)
b2 = bar(a2,'stacked'); grid on;
b2(1).Parent.Parent.Colormap = flipud(CM(4:7,:)); set(b2,'edgecolor','none');
ylim([0,100]); ylabel('Estimated incidence reduction in 2030 (%)');
legend('Coverage: PWID through OST and NSP', 'Coverage: all PWID','location','southeast');
set(gca,'XTicklabel',{'Current Ab\newlinetesting rates \newlinereplaced with \newlineRNA tests','Two-yearly \newlineRNA testing', 'Annual \newlineRNA testing', 'Six-monthly \newlineRNA testing'});
title(['\fontsize{14}2030 incidence reduction following treatment availability' char(10) '\it{}\fontsize{12}Antibody testing replaced with RNA testing among PWID'],'HorizontalAlignment','center');

a4 = a3; a4(2:end,2) = a3(2:end,2) + a3(2:end,1);
figure(3)
subplot(1,2,1)
b = bar(a4(2:end,2:end),'stacked'); grid on;
b(1).Parent.Parent.Colormap = flipud(CM1(3:8,:));
set(b,'edgecolor','none');
ylim([0,100]); ylabel('Estimated incidence reduction in 2030 (%)');
legend('50% RNA tested within 3 months', '60% ...', '70% ...', '80% ...',...
    '90% ...','100% RNA tested within 3 months', 'location','southeast');
set(gca,'XTicklabel',{'Two-yearly \newlineAb testing, \newlineall PWID',...
    'Annual \newlineAb testing, \newlineall PWID',...
    'Six-monthly \newlineAb testing, \newlineall PWID',...
    'Three-monthly \newlineAb testing, \newlineall PWID'});
title(['\it{}\fontsize{10}Changes to antibody testing and follow-up rates among PWID'],'HorizontalAlignment','center');
hold off;
subplot(1,2,2)
b2 = bar(a2(2:end,:),'stacked'); grid on;
b2(1).Parent.Parent.Colormap = flipud(CM1(3:8,:)); set(b2,'edgecolor','none');
ylim([0,100]); ylabel('Estimated incidence reduction in 2030 (%)');
legend('Coverage: PWID through OST and NSP', 'Coverage: all PWID','location','southeast');
set(gca,'XTicklabel',{'Two-yearly \newlineRNA testing', 'Annual \newlineRNA testing', 'Six-monthly \newlineRNA testing', 'Three-monthly \newlineRNA testing'});
title(['\it{}\fontsize{10}Antibody testing replaced with RNA testing among PWID'],'HorizontalAlignment','center');
axes('Position',[0 0 1 1],'Visible','off');
text(0.5,0.98,{'\fontsize{14}2030 incidence reduction following treatment availability'},'HorizontalAlignment','Center')



%% AUs case only

figure(4)
subplot(1,2,1)
b = bar(a3(1:end,:),'stacked'); grid on;
b(1).Parent.Parent.Colormap = flipud(CM1(3:8,:));
set(b,'edgecolor','none');
ylim([0,100]); ylabel('Estimated incidence reduction in 2030 (%)');
legend('Current average RNA follow-up rates in AUS (46%)',...
    '50% RNA tested within 3 months', '60% ...', '70% ...', '80% ...',...
    '90% ...','100% RNA tested within 3 months', 'location','southeast');
set(gca,'XTicklabel',{'Current \newlinetesting rates',...
    'Two-yearly \newlineAb testing, \newlineall PWID',...
    'Annual \newlineAb testing, \newlineall PWID',...
    'Six-monthly \newlineAb testing, \newlineall PWID',...
    'Three-monthly \newlineAb testing, \newlineall PWID'});
title(['\it{}\fontsize{10}Changes to antibody testing and follow-up rates among PWID'],'HorizontalAlignment','center');
hold off;
subplot(1,2,2)
b2 = bar(a2(1:end,:),'stacked'); grid on;
b2(1).Parent.Parent.Colormap = flipud(CM1(3:8,:)); set(b2,'edgecolor','none');
ylim([0,100]); ylabel('Estimated incidence reduction in 2030 (%)');
legend('Coverage: PWID through OST and NSP', 'Coverage: all PWID','location','southeast');
set(gca,'XTicklabel',{'Current Ab \newlinetesting rates \newlinereplaced with \newlineRNA tests',...
    'Two-yearly \newlineRNA testing', 'Annual \newlineRNA testing', 'Six-monthly \newlineRNA testing', 'Three-monthly \newlineRNA testing'});
title(['\it{}\fontsize{10}Antibody testing replaced with RNA testing among PWID'],'HorizontalAlignment','center');
axes('Position',[0 0 1 1],'Visible','off');
text(0.5,0.98,{'\fontsize{14}2030 incidence reduction following treatment availability'},'HorizontalAlignment','Center')




%% Calibrate to 30% prevalence
% scenario = 'empty';
% loaddata
% dt = 1/4;
% prev0 = [2015, 0.25;...
%     2014, 0.25;
%     2010, 0.25;
%     2008, 0.275;
%     2005, 0.3];
% alpha=alpha_old; % calibrate in pre-DAA era
% target_late=1; % target treatments to people with late-liver disease
% 
% [output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases,output_ost,output_nsp, output_diagnoses] = calibrate_optim_par(300, 20);
% save('C:\Users\Nick\Desktop\Matlab Sims\Testing\30prevcalibration1')
% 
% 
% % %% Calibrate to 60% prevalence
% scenario = 'empty';
% loaddata
% dt = 1/4;
% prev0 = [2015, 0.75;...
%     2014, 0.75;...
%     2010, 0.75;...
%     2008, 0.77;...
%     2005, 0.8];
% alpha=alpha_old; % calibrate in pre-DAA era
% target_late=1; % target treatments to people with late-liver disease
% 
% [output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases,output_ost,output_nsp, output_diagnoses] = calibrate_optim_par(300, 20);
% save('C:\Users\Nick\Desktop\Matlab Sims\Testing\60prevcalibration1')
