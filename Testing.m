clear all
global filename scenario sens target_late Tin Run start_year num_pops num_cascade num_age num_intervention num_engagement num_region dt ... %meta parameters
    mu mu_PWID mu_former exit_IDU r_relapse delta alpha alpha_old alpha_DAA  omega age_cohort P y0 t0...
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2... % disease progression parameters
    Q_sens q_svr q_svr_PWID q_treat  c_daa discount treat  ... % QALY and cost parameters
    p_complete total_PWID PWID0 age_mix prev0 r_inc_up followup...%PWID parameters
    prev0_MSM ... %MSM parameters
    incarceration_rate prison_length ... %Prison parameters
    cascade0 cascade0_PWID cascade0_MSM disease0 disease0_HIV cases0 ost0 nsp0 HCC0 diagnoses0 dem infect_factor infect_base progression_base imported... % data and calibration
    ost_duration  nsp_duration treat_projected target_inc target_death cascade_scale_time care RNAtesting... %intervention
    infect progression imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9 ost_enrollment nsp_enrollment... % calibtation parameters
    harm_reduction_coverage nsp_coverage ost_coverage prop_test diagnosed_risk_reduction

user=extractBetween(pwd,"Users\","\");
drive=extractBefore(pwd,":");
filename=strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\newanalysis");
filename2 = "G:\My Drive\Papers\02 Hepatitis C\17 Tanzania\Simulations\superPLHCV";

loaddata
%load('calibration_draftv2'); infect_base=infect; progression_base = progression;
dt = 1/4; % six-monthly time steps for burn-in / calibration perio 1950-2015
sens=1; %Number of runs in sensitivity analysis, sens=1 turns off feature
summary=zeros(6,12,sens);
%summary_HR = zeros(length(harm_reduction_range),length(summary(1,:,1)),sens);
followup = 1;
    

%% Sensitivity of parameters
for s=1:sens
    s
    scenario = 'empty';
    alpha=alpha_old; % calibrate in pre-DAA era
    target_late=1; % target treatments to people with late-liver disease
    infect_base = infect; progression_base = progression; treat = [0,0,0];
    %infect = infect_base; progression = progression_base;
    
    [output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases,output_ost,output_nsp, output_diagnoses] = ...
        calibrate_optim_par(200, 30);
    save(filename2); infect_base = infect; progression_base = progression; diagnosed_risk_reduction = 1; treat = [0,0,0];
    %load('calibration_test3'); infect_base=infect; progression_base = progression;
    
    
    
 %% Run model in prior to 2017
    [TT1,y1]=DE_track_age(Tin,y0,t0,treat);
    y1_end=reshape(y1(end,:,:,:,:,:,:,:), num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region,33);
    death2017 = (sum(sum(sum(sum(sum(sum(y1(end,:,:,:,:,:,:,22))))))) - sum(sum(sum(sum(sum(sum(y1(find(TT1>=66,1),:,:,:,:,:,:,22)))))))) / ...
        (TT1(end)-TT1(find(TT1>=66,1)));
    inc2017 = sum(sum(sum(sum(sum(sum(sum(y1(find(TT1>=66,1):end-1,:,:,:,:,:,:,27))))))));%/(TT1(find(TT1>=66,1)-1) - TT1(find(TT1>=65,1)));
    %y1_end(:,:,:,:,:,:,21:(27+6))=0;
    y1_end(:,6,:,:,:,:,1:20) = y1_end(:,6,:,:,:,:,1:20) + y1_end(:,8,:,:,:,:,1:20) + y1_end(:,10,:,:,:,:,1:20); 
    y1_end(:,8,:,:,:,:,1:20) = 0; y1_end(:,10,:,:,:,:,1:20) = 0; % moving failed to be re-eligible with DAAs
    
    sum(sum(sum(sum(sum(y1(find(TT1>=Tin-(2016-prev0(1,1)),1),1,:,:,:,:,1,6:20))))))./sum(sum(sum(sum(sum(y1(find(TT1>=Tin-(2016-prev0(1,1)),1),1,:,:,:,:,1,1:20))))))
    
    if sens > 1
        perturb_parameters
    end
    %% Baseline: Current standard of care with no scaled up treatment
    scenario = 'base'; %Current level of community care
    alpha = alpha_DAA;
    prop_test = 1;
     
    [TT2,y2]=DE_track_age(Run,y1_end,TT1,treat);
    [ycomb_noage, summary(1,:,s), tr, tr_] = gather_outputs(y1,y2,TT2);
    
 
    %%  Scenario 1: Scaled harm reduction
    scenario = 'current'; cal = 0;
    harm_reduction_range = [0,-1,0.1:0.1:0.5];
    
    for h = 1:length(harm_reduction_range)
        if harm_reduction_range(h)<0
            nsp_coverage = 0.01;
            ost_coverage = 0.1;
        else
            nsp_coverage = harm_reduction_range(h);
            ost_coverage = harm_reduction_range(h);
        end
        [TT2_treat,y2_treat]=DE_track_age(Run,y1_end,TT1,treat);
        [ycomb2_noage, summary(2,:,s), tr2, tr2_] = gather_outputs(y1,y2_treat,TT2_treat);
        
        TT_ = TT2_treat(TT2_treat-(Tin-10)>0)-(Tin-10);
        y2comb=[y1(1:end-1,:,:,:,:,:,:,:);y2_treat];
        y2_noage=reshape(sum(sum(sum(sum(y2_treat,4),5),6),7),size(y2_treat,1),num_pops, num_cascade,27+6); %Reshape to sum over age stratification, interventions, engagement and region
        
        squeeze(sum(sum(sum(sum(sum(sum(y2_treat(:,1,:,:,:,:,:,1:20),2),3),4),6),7),8)); %Reshape to sum over age stratification, interventions, engagement and region
        
        %         for i = Tin-10:Tin
        %             prev_HR(i,h) = 100*sum(sum(sum(sum(sum(sum(y1(find(TT1>=i,1),1,:,:,:,:,:,[6,12:20])))))))./...
        %                 sum(sum(sum(sum(sum(sum(y1(find(TT1>=i,1),1,:,:,:,:,:,1:20)))))));
        %             inc_HR(i,h) = sum(sum(sum(sum(sum(sum(sum(y1(find(TT1>=i-1,1):find(TT1>=i,1)-1,:,:,:,:,:,:,27))))))));
        %         end
        
        for i = (Tin-10):(Tin+Run)
            inc_HR(i-Tin+11,h,s) = (sum(sum(sum(ycomb2_noage(find(TT2_treat>=i-1,1):find(TT2_treat>=i,1)-1,:,:,27)))));
            prev_HR(i-Tin+11,h,s) = 100*sum(sum(ycomb2_noage(find(TT2_treat>=i,1),1,:,[6,12:20])))./...
                sum(sum(ycomb2_noage(find(TT2_treat>=i,1),1,:,1:20)));
            death_HR(i-Tin+11,h,s) = (sum(sum(ycomb2_noage(find(TT2_treat>=i,1)-1,:,:,22))) - sum(sum(ycomb2_noage(find(TT2_treat>=i-1,1),:,:,22))))...
                ./(TT2_treat(find(TT2_treat>=i,1))-TT2_treat(find(TT2_treat>=i-1,1)));
            [costs_HR(i-Tin+11,h,s),~,~,~,c_HRonly(i-Tin+11,h,s)] = Costs_age(...
                ycomb2_noage(find(TT2_treat>=i-1,1):find(TT2_treat>=i,1),:,:,:,:,:,:,:),...
                TT2_treat(find(TT2_treat>=i-1,1):find(TT2_treat>=i,1))-Tin,...
                y2comb(find(TT2_treat>=i-1,1):find(TT2_treat>=i,1),:,:,:,:,:,:,:));
        end
        summary_HR(h,:,s) = summary(2,:,s);
    end
    nsp_coverage = 0.01;
    ost_coverage = 0.1;
    
    %%  Scenario 2: Ab + RNA (standard testing) to reach 90% diagnosed
    scenario = 'current';
    target_late=1; % Target PWID
    alpha = alpha_DAA;
    progression = progression_base;
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    progression(1,5,2,1) = 1/(30/365); progression(2,5,2,1) = 1/(30/365); progression(3,5,2,1) = 1/(30/365); 
    progression(1,6,2,1) = 1/(30/365); progression(2,6,2,1) = 1/(30/365); progression(3,6,2,1) = 1/(30/365); 
    range_test = [0, 0.5]*2; % divided by 2 to assume that infection happens midway between tests
    range_followup = [0.74];
    range_treat = [0,-0.8]; % 80% of diagnosed treated according to WHO
    %range_test = [4]*2; % divided by 2 to assume that infection happens midway between tests
    range_diagnosed_risk_reduction = [0:0.1:0.1]; % risk reduction when diagnosed
    prop_test_range = [0.7]; % 80% coverage
    harm_reduction_range = [0,-1,0.1:0.1:0.5];
    for i = 1:length(range_test)
        for j = 1:length(range_treat)
            for k =1:length(harm_reduction_range)
                prop_test = prop_test_range(1);
                if harm_reduction_range(k)<0
                    nsp_coverage = 0.01;
                    ost_coverage = 0.1;
                else
                    nsp_coverage = harm_reduction_range(k);
                    ost_coverage = harm_reduction_range(k);
                end
                if range_test(i) > 0
                    progression(1,1,2,1) = range_test(i);
                    progression(2,1,2,1) = range_test(i);
                    %progression(3,1,2,1) = range_test(i);
                else
                    progression(1,1,2,1) = progression_base(1,1,2,1);
                    progression(2,1,2,1) = progression_base(2,1,2,1);
                    progression(3,1,2,1) = progression_base(3,1,2,1);
                    followup = 1;
                end
                if range_followup(1) > 0
                    followup = range_followup(1);%/0.46;
                    progression(2,2,2,1) = 1/0.25; progression(1,2,2,1) = 1 / 0.25; progression(3,2,2,1) = 1 / 0.25;
                else
                    progression(1,2,2,1) = progression_base(1,2,2,1);
                    progression(2,2,2,1) = progression_base(2,2,2,1);
                    progression(3,2,2,1) = progression_base(3,2,2,1);
                    followup = 1;
                end
                treat_no = [range_treat(j),0,0];
                y1_end_sim = y1_end;
                y1_end_sim(:,:,:,:,1,:,:) = (1-prop_test)* y1_end(:,:,:,:,2,:,:);
                y1_end_sim(:,:,:,:,2,:,:) = prop_test * y1_end(:,:,:,:,2,:,:);
                y1_end_sim(:,1,:,:,1,:,1:20) = y1_end_sim(:,1,:,:,1,:,1:20) + sum(y1_end_sim(:,2:end,:,:,1,:,1:20),2);
                y1_end_sim(:,2:end,:,:,1,:,1:20) = 0;
                [TT2_treat5,y2_treat5]=DE_track_age(Run,y1_end_sim,TT1,treat_no);
                [ycomb5_noage, summary(5,:,s), tr5, tr5_] = gather_outputs(y1,y2_treat5,TT2_treat5);
                
                for t = (Tin-10):(Tin+Run)
                    inc_test(t-Tin+11,i,j,k) = (sum(sum(sum(ycomb5_noage(find(TT2_treat5>=t-1,1):find(TT2_treat5>=t,1)-1,:,:,27)))));
                    prev_test(t-Tin+11,i,j,k) = 100*sum(sum(ycomb5_noage(find(TT2_treat5>=t,1),1,:,[6,12:20])))./...
                        sum(sum(ycomb5_noage(find(TT2_treat5>=t,1),1,:,1:20)));
%                     diag_test(t-Tin+11,i,j,k) = 100*sum(sum(sum(ycomb5_noage(find(TT2_treat5>=t,1),:,3:end,[6,12:20]))))./...
%                         sum(sum(sum(ycomb5_noage(find(TT2_treat5>=t,1),:,:,[6,12:20]))));
                    diag_test(t-Tin+11,i,j,k) = 100*(1 - sum(sum(sum(ycomb5_noage(find(TT2_treat5>=t,1),:,1,[6,12:20]))))./...
                        sum(sum(sum(ycomb5_noage(find(TT2_treat5>=Tin,1),:,:,[6,12:20]))))); % 1-undiagnosed at time t divided by total in 2015
                    death_test(t-Tin+11,i,j,k) = (sum(sum(ycomb5_noage(find(TT2_treat5>=t,1)-1,:,:,22))) - sum(sum(ycomb5_noage(find(TT2_treat5>=t-1,1),:,:,22))))...
                        ./(TT2_treat5(find(TT2_treat5>=t,1))-TT2_treat5(find(TT2_treat5>=t-1,1)));
                end
                summary_test(i,j,k,:,s) = summary(5,:,s);
            end
        end
    end
    
    %%  Scenario 3:  DBS HCVcAg testing
    scenario = 'DBS_HCVcAg';
    target_late=0; % Target PWID
    alpha = alpha_DAA;
    progression = progression_base;
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    progression(1,2,2,1) = 50; progression(2,2,2,1) = 50; progression(3,2,2,1) = 50; % perfect "RNA" followup
    progression(1,5,2,1) = 1/(30/365); progression(2,5,2,1) = 1/(30/365); progression(3,5,2,1) = 1/(30/365); 
    progression(1,6,2,1) = 1/(30/365); progression(2,6,2,1) = 1/(30/365); progression(3,6,2,1) = 1/(30/365); 
    range_followup = [0.76]; % Reduced test sensitivity
    range_treat = [0,-0.8];
    range_test = [0,0.5]*2; % divided by 2 to assume that infection happens midway between tests
    range_diagnosed_risk_reduction = [0:0.1:0.1]; % risk reduction when diagnosed
    prop_test_range = [0.9];
    harm_reduction_range = [0,-1,0.1:0.1:0.5];
    for i = 1:length(range_test)
        for j = 1:length(range_treat)
            for k =1:length(harm_reduction_range)
                prop_test = prop_test_range(1);
                followup = range_followup(1);
                if harm_reduction_range(k)<0
                    nsp_coverage = 0.01;
                    ost_coverage = 0.1;
                else
                    nsp_coverage = harm_reduction_range(k);
                    ost_coverage = harm_reduction_range(k);
                end
                if range_test(i) > 0
                    progression(1,1,2,1) = range_test(i);
                    progression(2,1,2,1) = range_test(i);
                    %progression(3,1,2,1) = range_test(i);
                else
                    progression(1,1,2,1) = progression_base(1,1,2,1);
                    progression(2,1,2,1) = progression_base(2,1,2,1);
                    progression(3,1,2,1) = progression_base(3,1,2,1);
                    followup = 1;
                end
                treat_no = [range_treat(j),0,0];
                y1_end_sim = y1_end;
                y1_end_sim(:,:,:,:,1,:,:) = (1-prop_test)* y1_end(:,:,:,:,2,:,:);
                y1_end_sim(:,:,:,:,2,:,:) = prop_test * y1_end(:,:,:,:,2,:,:);
                y1_end_sim(:,1,:,:,1,:,1:20) = y1_end_sim(:,1,:,:,1,:,1:20) + sum(y1_end_sim(:,2:end,:,:,1,:,1:20),2);
                y1_end_sim(:,2:end,:,:,1,:,1:20) = 0;
                [TT2_treat3,y2_treat3]=DE_track_age(Run,y1_end_sim,TT1,treat_no);
                [ycomb3_noage, summary(3,:,s), tr3, tr3_] = gather_outputs(y1,y2_treat3,TT2_treat3);
                
                for t = (Tin-10):(Tin+Run)
                    inc_DBS(t-Tin+11,i,j,k) = (sum(sum(sum(ycomb3_noage(find(TT2_treat3>=t-1,1):find(TT2_treat3>=t,1)-1,:,:,27)))));
                    prev_DBS(t-Tin+11,i,j,k) = 100*sum(sum(ycomb3_noage(find(TT2_treat3>=t,1),1,:,[6,12:20])))./...
                        sum(sum(ycomb3_noage(find(TT2_treat3>=t,1),1,:,1:20)));
%                     diag_DBS(t-Tin+11,i,j,k) = 100*sum(sum(sum(ycomb3_noage(find(TT2_treat3>=t,1),:,2:end,[6,12:20]))))./...
%                         sum(sum(sum(ycomb3_noage(find(TT2_treat3>=t,1),:,:,[6,12:20]))));
                    diag_DBS(t-Tin+11,i,j,k) = 100*(1 - sum(sum(sum(ycomb3_noage(find(TT2_treat3>=t,1),:,1,[6,12:20]))))./...
                        sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin,1),:,:,[6,12:20]))))); % 1-undiagnosed at time t divided by total in 2015
                    death_DBS(t-Tin+11,i,j,k) = (sum(sum(ycomb3_noage(find(TT2_treat3>=t,1)-1,:,:,22))) - sum(sum(ycomb3_noage(find(TT2_treat3>=t-1,1),:,:,22))))...
                        ./(TT2_treat3(find(TT2_treat3>=t,1))-TT2_treat3(find(TT2_treat3>=t-1,1)));
                end
                summary_DBS(i,j,k,:,s) = summary(3,:,s);
            end
        end
    end
    
    
    %%  Scenario 4: serum HCVcAg testing
    scenario = 'serum_HCVcAg';
    target_late=0; % Target PWID
    alpha = alpha_DAA;
    progression = progression_base;
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    progression(1,2,2,1) = 50; progression(2,2,2,1) = 50; progression(3,2,2,1) = 50; % perfect "RNA" followup
    progression(1,5,2,1) = 1/(30/365); progression(2,5,2,1) = 1/(30/365); progression(3,5,2,1) = 1/(30/365); 
    progression(1,6,2,1) = 1/(30/365); progression(2,6,2,1) = 1/(30/365); progression(3,6,2,1) = 1/(30/365); 
    range_followup = [1];
    range_treat = [0,-0.8];
    range_test = [0,0.5]*2; % divided by 2 to assume that infection happens midway between tests
    range_diagnosed_risk_reduction = [0:0.1:0.1]; % risk reduction when diagnosed
    prop_test_range = [0.7];
    harm_reduction_range = [0,-1,0.1:0.1:0.5];
    for i = 1:length(range_test)
        for j = 1:length(range_treat)
            for k =1:length(harm_reduction_range)
                prop_test = prop_test_range(1);
                followup = range_followup(1);
                if harm_reduction_range(k)<0
                    nsp_coverage = 0.01;
                    ost_coverage = 0.1;
                else
                    nsp_coverage = harm_reduction_range(k);
                    ost_coverage = harm_reduction_range(k);
                end
                if range_test(i) > 0
                    progression(1,1,2,1) = range_test(i);
                    progression(2,1,2,1) = range_test(i);
                    %progression(3,1,2,1) = range_test(i);
                else
                    progression(1,1,2,1) = progression_base(1,1,2,1);
                    progression(2,1,2,1) = progression_base(2,1,2,1);
                    progression(3,1,2,1) = progression_base(3,1,2,1);
                    followup = 1;
                end
                treat = [range_treat(j),0,0];
                y1_end_sim = y1_end;
                y1_end_sim(:,:,:,:,1,:,:) = (1-prop_test)* y1_end(:,:,:,:,2,:,:);
                y1_end_sim(:,:,:,:,2,:,:) = prop_test * y1_end(:,:,:,:,2,:,:);
                y1_end_sim(:,1,:,:,1,:,1:20) = y1_end_sim(:,1,:,:,1,:,1:20) + sum(y1_end_sim(:,2:end,:,:,1,:,1:20),2);
                y1_end_sim(:,2:end,:,:,1,:,1:20) = 0;
                [TT2_treat6,y2_treat6]=DE_track_age(Run,y1_end_sim,TT1,treat);
                [ycomb6_noage, summary(6,:,s), tr6, tr6_] = gather_outputs(y1,y2_treat6,TT2_treat6);
                
                for t = (Tin-10):(Tin+Run)
                    inc_serum(t-Tin+11,i,j,k) = (sum(sum(sum(ycomb6_noage(find(TT2_treat6>=t-1,1):find(TT2_treat6>=t,1)-1,:,:,27)))));
                    prev_serum(t-Tin+11,i,j,k) = 100*sum(sum(ycomb6_noage(find(TT2_treat6>=t,1),1,:,[6,12:20])))./...
                        sum(sum(ycomb6_noage(find(TT2_treat6>=t,1),1,:,1:20)));
%                     diag_serum(t-Tin+11,i,j,k) = 100*sum(sum(sum(ycomb6_noage(find(TT2_treat6>=t,1),:,2:end,[6,12:20]))))./...
%                         sum(sum(sum(ycomb6_noage(find(TT2_treat6>=t,1),:,:,[6,12:20]))));
                    diag_serum(t-Tin+11,i,j,k) = 100*(1 - sum(sum(sum(ycomb6_noage(find(TT2_treat6>=t,1),:,1,[6,12:20]))))./...
                        sum(sum(sum(ycomb6_noage(find(TT2_treat6>=Tin,1),:,:,[6,12:20]))))); % 1-undiagnosed at time t divided by total in 2015
                    death_serum(t-Tin+11,i,j,k) = (sum(sum(ycomb6_noage(find(TT2_treat6>=t,1)-1,:,:,22))) - sum(sum(ycomb6_noage(find(TT2_treat6>=t-1,1),:,:,22))))...
                        ./(TT2_treat6(find(TT2_treat6>=t,1))-TT2_treat6(find(TT2_treat6>=t-1,1)));
                end
                summary_serum(i,j,k,:,s) = summary(6,:,s);
            end
        end
    end
    
    
    %%  Scenario 5: Diagnosed risk reduction
    scenario = 'DBS_HCVcAg';
    target_late=0; % Target PWID
    alpha = alpha_DAA;
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    progression(1,2,2,1) = 50; progression(2,2,2,1) = 50; progression(3,2,2,1) = 50; % perfect "RNA" followup
    range_followup = [0.74];
    range_test = [0,0.5]; % divided by 2 to assume that infection happens midway between tests
    range_diagnosed_risk_reduction = [0:0.1:0.1]; % risk reduction when diagnosed
    prop_test_range = [0.9];
    followup = 0.76; %DBS sensitivity
    for i = 1:length(range_test)
        for j = 1:length(range_diagnosed_risk_reduction)
            for k =1:length(prop_test_range)
                prop_test = prop_test_range(k);
                diagnosed_risk_reduction = (1-range_diagnosed_risk_reduction(j));
                if range_test(i) > 0
                    progression(1,1,2,1) = range_test(i);
                    progression(2,1,2,1) = range_test(i);
                    progression(3,1,2,1) = range_test(i);
                else
                    progression(1,1,2,1) = progression_base(1,1,2,1);
                    progression(2,1,2,1) = progression_base(2,1,2,1);
                    progression(3,1,2,1) = progression_base(3,1,2,1);
                    followup = 1;
                end
                y1_end_sim = y1_end;
                y1_end_sim(:,:,:,:,1,:,:) = (1-prop_test)* y1_end(:,:,:,:,2,:,:);
                y1_end_sim(:,:,:,:,2,:,:) = prop_test * y1_end(:,:,:,:,2,:,:);
                y1_end_sim(:,1,:,:,1,:,1:20) = y1_end_sim(:,1,:,:,1,:,1:20) + sum(y1_end_sim(:,2:end,:,:,1,:,1:20),2);
                y1_end_sim(:,2:end,:,:,1,:,1:20) = 0;
                [TT2_treat4,y2_treat4]=DE_track_age(Run,y1_end_sim,TT1,treat);
                [ycomb4_noage, summary(4,:,s), tr4, tr4_] = gather_outputs(y1,y2_treat4,TT2_treat4);
                
                for t = (Tin-10):(Tin+Run)
                    inc_DBS_rr(t-Tin+11,i,j,k) = (sum(sum(sum(ycomb4_noage(find(TT2_treat4>=t-1,1):find(TT2_treat4>=t,1)-1,:,:,27)))));
                    prev_DBS_rr(t-Tin+11,i,j,k) = 100*sum(sum(ycomb4_noage(find(TT2_treat4>=t,1),1,:,[6,12:20])))./...
                        sum(sum(ycomb4_noage(find(TT2_treat4>=t,1),1,:,1:20)));
                    diag_DBS_rr(t-Tin+11,i,j,k) = 100*sum(sum(sum(ycomb4_noage(find(TT2_treat4>=t,1),:,2:end,[6,12:20]))))./...
                        sum(sum(sum(ycomb4_noage(find(TT2_treat4>=t,1),:,:,[6,12:20]))));
                end
                summary_DBS_rr(i,j,k,:,s) = summary(4,:,s);
            end
        end
    end

    
   
    %% Collect outputs
    years=60:1:81;
    for t=2:length(years)
        inc_year_sens(:,t-1,s)=[sum(sum(sum(ycomb_noage(find(TT2>=years(t-1),1):find(TT2>=years(t),1)-1,:,:,27))));...
            sum(sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1):find(TT2_treat>=years(t),1)-1,:,:,27))));...
            sum(sum(sum(ycomb3_noage(find(TT2_treat3>=years(t-1),1):find(TT2_treat3>=years(t),1)-1,:,:,27))));...
            sum(sum(sum(ycomb4_noage(find(TT2_treat4>=years(t-1),1):find(TT2_treat4>=years(t),1)-1,:,:,27))));...
            sum(sum(sum(ycomb5_noage(find(TT2_treat5>=years(t-1),1):find(TT2_treat5>=years(t),1)-1,:,:,27))));...
            sum(sum(sum(ycomb6_noage(find(TT2_treat6>=years(t-1),1):find(TT2_treat6>=years(t),1)-1,:,:,27))))];%...
            %./[TT2(find(TT2>=years(t),1)-1)-TT2(find(TT2>=years(t-1),1));TT2_treat(find(TT2_treat>=years(t),1)-1)-TT2_treat(find(TT2_treat>=years(t-1),1));...
            %TT2_treat3(find(TT2_treat3>=years(t),1)-1)-TT2_treat3(find(TT2_treat3>=years(t-1),1));TT2_treat4(find(TT2_treat4>=years(t),1)-1)-TT2_treat4(find(TT2_treat4>=years(t-1),1));...
            %TT2_treat5(find(TT2_treat5>=years(t),1)-1)-TT2_treat5(find(TT2_treat5>=years(t-1),1));TT2_treat6(find(TT2_treat6>=years(t),1)-1)-TT2_treat6(find(TT2_treat6>=years(t-1),1))];
        
        death_year_sens(:,t-1,s)=[sum(sum(ycomb_noage(find(TT2>=years(t),1)-1,:,:,22))) - sum(sum(ycomb_noage(find(TT2>=years(t-1),1),:,:,22)));...
            sum(sum(ycomb2_noage(find(TT2_treat>=years(t),1)-1,:,:,22))) - sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1),:,:,22)));...
            sum(sum(ycomb3_noage(find(TT2_treat3>=years(t),1)-1,:,:,22))) - sum(sum(ycomb3_noage(find(TT2_treat3>=years(t-1),1),:,:,22)));...
            sum(sum(ycomb4_noage(find(TT2_treat4>=years(t),1)-1,:,:,22))) - sum(sum(ycomb4_noage(find(TT2_treat4>=years(t-1),1),:,:,22)));...
            sum(sum(ycomb5_noage(find(TT2_treat5>=years(t),1)-1,:,:,22))) - sum(sum(ycomb5_noage(find(TT2_treat5>=years(t-1),1),:,:,22)));...
            sum(sum(ycomb6_noage(find(TT2_treat6>=years(t),1)-1,:,:,22))) - sum(sum(ycomb6_noage(find(TT2_treat6>=years(t-1),1),:,:,22)))]...
            ./...
            [TT2(find(TT2>=years(t),1))-TT2(find(TT2>=years(t-1),1));...
            TT2_treat(find(TT2_treat>=years(t),1)-1)-TT2_treat(find(TT2_treat>=years(t-1),1));...
            TT2_treat3(find(TT2_treat3>=years(t),1)-1)-TT2_treat3(find(TT2_treat3>=years(t-1),1));...
            TT2_treat4(find(TT2_treat4>=years(t),1)-1)-TT2_treat4(find(TT2_treat4>=years(t-1),1));...
            TT2_treat5(find(TT2_treat5>=years(t),1)-1)-TT2_treat5(find(TT2_treat5>=years(t-1),1));...
            TT2_treat6(find(TT2_treat6>=years(t),1)-1)-TT2_treat6(find(TT2_treat6>=years(t-1),1))];
        
        treatment_alloc_sens(:,t-1,s)=[sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1):find(TT2_treat>years(t),1)-1,1,:,23)));...% PWID early liver disease treatments
            sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1):find(TT2_treat>years(t),1)-1,1,:,25)));... PWID late liver disease treatments
            sum(sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1):find(TT2_treat>years(t),1)-1,2:3,:,23))));... %Former PWID early liver disease stage
            sum(sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1):find(TT2_treat>years(t),1)-1,2:3,:,25))))];%... %Former PWID late liver disease stage
        %./repmat(TT2_treatlate(find(TT2_treatlate>=years(t),1)-1)-TT2_treatlate(find(TT2_treatlate>=years(t-1),1)),4,1);
    end
    years2=50:1:81;
    for t=2:length(years2)
        HCC_year_sens(:,t-1,s)=[sum(sum(sum(ycomb_noage(find(TT2>=years2(t-1),1):find(TT2>=years2(t),1)-1,:,:,24))));...
            sum(sum(sum(ycomb2_noage(find(TT2_treat>=years2(t-1),1):find(TT2_treat>=years2(t),1)-1,:,:,24))));...
            sum(sum(sum(ycomb3_noage(find(TT2_treat3>=years2(t-1),1):find(TT2_treat3>=years2(t),1)-1,:,:,24))));...
            sum(sum(sum(ycomb4_noage(find(TT2_treat4>=years2(t-1),1):find(TT2_treat4>=years2(t),1)-1,:,:,24))));...
            sum(sum(sum(ycomb5_noage(find(TT2_treat5>=years2(t-1),1):find(TT2_treat5>=years2(t),1)-1,:,:,24))));...
            sum(sum(sum(ycomb6_noage(find(TT2_treat6>=years2(t-1),1):find(TT2_treat6>=years2(t),1)-1,:,:,24))))]...
            ./[TT2(find(TT2>=years2(t),1)-1)-TT2(find(TT2>=years2(t-1),1));TT2_treat(find(TT2_treat>=years2(t),1)-1)-TT2_treat(find(TT2_treat>=years2(t-1),1));...
            TT2_treat3(find(TT2_treat3>=years2(t),1)-1)-TT2_treat3(find(TT2_treat3>=years2(t-1),1));TT2_treat4(find(TT2_treat4>=years2(t),1)-1)-TT2_treat4(find(TT2_treat4>=years2(t-1),1));...
            TT2_treat5(find(TT2_treat5>=years2(t),1)-1)-TT2_treat5(find(TT2_treat5>=years2(t-1),1));TT2_treat6(find(TT2_treat6>=years2(t),1)-1)-TT2_treat6(find(TT2_treat6>=years2(t-1),1))];
        
        diagnosed_year_sens(:,t-1,s)=[sum(sum(sum(ycomb_noage(find(TT2>=years2(t-1),1):find(TT2>=years2(t),1)-1,:,:,29))));...
            sum(sum(sum(ycomb2_noage(find(TT2_treat>=years2(t-1),1):find(TT2_treat>=years2(t),1)-1,:,:,29))));...
            sum(sum(sum(ycomb3_noage(find(TT2_treat3>=years2(t-1),1):find(TT2_treat3>=years2(t),1)-1,:,:,29))));...
            sum(sum(sum(ycomb4_noage(find(TT2_treat4>=years2(t-1),1):find(TT2_treat4>=years2(t),1)-1,:,:,29))));...
            sum(sum(sum(ycomb5_noage(find(TT2_treat5>=years2(t-1),1):find(TT2_treat5>=years2(t),1)-1,:,:,29))));...
            sum(sum(sum(ycomb6_noage(find(TT2_treat6>=years2(t-1),1):find(TT2_treat6>=years2(t),1)-1,:,:,29))))]...
            ./[TT2(find(TT2>=years2(t),1)-1)-TT2(find(TT2>=years2(t-1),1));TT2_treat(find(TT2_treat>=years2(t),1))-TT2_treat(find(TT2_treat>=years2(t-1),1));...
            TT2_treat3(find(TT2_treat3>=years2(t),1)-1)-TT2_treat3(find(TT2_treat3>=years2(t-1),1));TT2_treat4(find(TT2_treat4>=years2(t),1)-1)-TT2_treat4(find(TT2_treat4>=years2(t-1),1));...
            TT2_treat5(find(TT2_treat5>=years2(t),1)-1)-TT2_treat5(find(TT2_treat5>=years2(t-1),1));TT2_treat6(find(TT2_treat6>=years2(t),1)-1)-TT2_treat6(find(TT2_treat6>=years2(t-1),1))];
    end
    for t = 1:length(years)
        Total_HCV_sens(t,:,s) = [sum(sum(sum(ycomb_noage(find(TT2>=years(t),1),:,:,6:20)))),sum(sum(sum(ycomb2_noage(find(TT2_treat>=years(t),1),:,:,6:20)))),...
            sum(sum(sum(ycomb5_noage(find(TT2_treat5>=years(t),1),:,:,6:20)))),sum(sum(sum(ycomb6_noage(find(TT2_treat6>=years(t),1),:,:,6:20)))),...
            sum(sum(sum(ycomb3_noage(find(TT2_treat3>=years(t),1),:,:,6:20)))),sum(sum(sum(ycomb4_noage(find(TT2_treat4>=years(t),1),:,:,6:20))))];
    end
    Total_HCV_sens(:,:,s) = round(Total_HCV_sens(:,:,s) / 1000);
    
    %inc_year_sens(:,1,s) = [inc2015;inc2015;inc2015;inc2015;inc2015;inc2015];
    death_year_sens(:,6,s) = [death2017;death2017;death2017;death2017;death2017;death2017];
    total_treat(:,:,s) = [tr; tr2; tr3; tr4; tr5; tr6];
    total_treat_2030(:,:,s) = [tr_; tr2_; tr3_; tr4_; tr5_; tr6_];
    
end


inc_year=mean(inc_year_sens(:,:,:),3);%./repmat(inc_year_sens(1,1,:),6,length(years)-1),3);
%inc_year(5,:)=mean(inc_year_sens(5,:,treat_scaleup(1,:)>=0)./repmat(inc_year_sens(1,1,treat_scaleup(1,:)>=0),1,length(years)-1),3);

%inc_year = mean(inc_year_sens(:,:,:),3);
inc_year_std=prctile(inc_year_sens(:,:,:),[5,95],3);%./repmat(inc_year_sens(1,1,:),6,length(years)-1)),[5,95],3);
%inc_year_std(5,:,:)=prctile((inc_year_sens(5,:,treat_scaleup(1,:)>=0)./repmat(inc_year_sens(1,1,treat_scaleup(1,:)>=0),1,length(years)-1)),[5,95],3);
deaths15_30=mean(sum(death_year_sens(:,[1:15],:),2),3);
deaths15_30_std=prctile(sum(death_year_sens(:,[1:15],:),2),[5,95],3);
death_year=mean(death_year_sens(:,:,:),3);
death_year_std=prctile((death_year_sens(:,:,:)./repmat(death_year_sens(1,1,:),6,length(years)-1)),[5,95],3);
treatment_alloc=mean(treatment_alloc_sens,3);
treatment_alloc_std=std(treatment_alloc_sens,0,3);
HCC_year=mean(HCC_year_sens(:,:,:),3);
diagnosed_year=mean(diagnosed_year_sens(:,:,:),3);
% treat_scaleup_summary(:,1) = mean(treat_scaleup,2); treat_scaleup_summary(1,1) = mean(treat_scaleup(1,treat_scaleup(1,:)>=0),2);
% treat_scaleup_summary(:,2:3) = prctile(treat_scaleup,[25,75],2); treat_scaleup_summary(1,2:3) = prctile(treat_scaleup(1,treat_scaleup(1,:)>=0),[25,75],2);
% R0_summary(1,:) = mean(R0_scens,1); R0_summary(2:3,:) = prctile(R0_scens, [5,95],1);

for i = 1:length(Total_HCV_sens(:,1,1))
    for j = 1:6
%         if i==5
%             Total_HCV(i,j,:) = prctile(Total_HCV_sens(i,j,treat_scaleup(1,:)>=0),[5,95]);
%         else
            Total_HCV(i,j,:) = prctile(Total_HCV_sens(i,j,:),[5,95]);
%         end
    end
end

summary(1:6,5,:)= 100*(summary(1:6,5,:)- repmat(inc_year_sens(1,1,:),[6,1]))./repmat(inc_year_sens(1,1,:),[6,1]); %Percent reduction in incidence
summary(1:6,6,:)= 100*(summary(1:6,6,:)- repmat(death_year_sens(1,1,:),[6,1]))./repmat(death_year_sens(1,1,:),[6,1]); %Percent reduction in mortality

%Additional costs/QALYs compared to baseline
summary2=zeros(6,12,sens);
summary2(1,:,:)=summary(1,:,:);
summary2(2:6,:,:)=summary(2:6,:,:)-repmat(summary(1,:,:),[5,1,1]);
summary2(2:6,4,:)=summary(2:6,4,:);

% Confidence intervals on outputs from uncertainty analysis
margins=zeros(6,12,3); 
margins(:,:,1)=mean(summary(:,:,:),3); %margins(5,:,1)=mean(summary(5,:,treat_scaleup(1,:)>=0),3);
margins2=zeros(6,12,3); 
margins2(:,:,1)=mean(summary2(:,:,:),3); %margins2(5,:,1)=mean(summary2(5,:,treat_scaleup(1,:)>=0),3);
total_treat_summary = zeros(6,4,3); 
total_treat_summary(:,:,1) = mean(total_treat(:,:,:),3); %total_treat_summary(5,:,1) = mean(total_treat(5,:,treat_scaleup(1,:)>=0),3);
total_treat_2030_summary = zeros(6,4,3); 
total_treat_2030_summary(:,:,1) = mean(total_treat_2030(:,:,:),3); %total_treat_2030_summary(5,:,1) = mean(total_treat_2030(5,:,treat_scaleup(1,:)>=0),3);

for i=1:6
    for j=1:12
%         if i==5
%             margins(i,j,2)=prctile(summary(i,j,treat_scaleup(1,:)>=0),5); %Lower bound
%             margins(i,j,3)=prctile(summary(i,j,treat_scaleup(1,:)>=0),95); %Upper bound
%         
%             margins2(i,j,2)=prctile(summary2(i,j,treat_scaleup(1,:)>=0),5); %Lower bound
%             margins2(i,j,3)=prctile(summary2(i,j,treat_scaleup(1,:)>=0),95); %Upper bound
%         else
            margins(i,j,2)=prctile(summary(i,j,:),5); %Lower bound
            margins(i,j,3)=prctile(summary(i,j,:),95); %Upper bound
        
            margins2(i,j,2)=prctile(summary2(i,j,:),5); %Lower bound
            margins2(i,j,3)=prctile(summary2(i,j,:),95); %Upper bound
%         end
    end
    for j=1:4
%         if i==5
%             total_treat_summary(i,j,2)=prctile(total_treat(i,j,treat_scaleup(1,:)>=0),25); %Lower bound
%             total_treat_summary(i,j,3)=prctile(total_treat(i,j,treat_scaleup(1,:)>=0),75); %Upper bound
%             total_treat_2030_summary(i,j,2)=prctile(total_treat_2030(i,j,treat_scaleup(1,:)>=0),25); %Lower bound
%             total_treat_2030_summary(i,j,3)=prctile(total_treat_2030(i,j,treat_scaleup(1,:)>=0),75); %Upper bound
%         else
            total_treat_summary(i,j,2)=prctile(total_treat(i,j,:),25); %Lower bound
            total_treat_summary(i,j,3)=prctile(total_treat(i,j,:),75); %Upper bound
            total_treat_2030_summary(i,j,2)=prctile(total_treat_2030(i,j,:),25); %Lower bound
            total_treat_2030_summary(i,j,3)=prctile(total_treat_2030(i,j,:),75); %Upper bound
%         end
    end
end
save(filename)
charts_Tanz
% %Table for paper
% margins_point = load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\baseline_point','margins');
% total_treat_point = load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\baseline_point','total_treat_summary');
% total_treat_2030_summary_point = load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\baseline_point','total_treat_2030_summary');
% treat_scaleup_point = load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\baseline_point','treat_scaleup');
% inc_year = load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\point','inc_year');
% R0_summary_point = load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\baseline_point','R0_summary');
% margins_point = struct2array(margins_point); total_treat_point = struct2array(total_treat_point);
% treat_scaleup_point = struct2array(treat_scaleup_point);
% total_treat_2030_summary_point = struct2array(total_treat_2030_summary_point);
% inc_year = struct2array(inc_year); R0_summary_point = struct2array(R0_summary_point);

summary_HR_point = struct2array(load('Base_draftv2','summary_HR'));
inc_HR_point = struct2array(load('Base_draftv2','inc_HR'));
death_year_sens_point = struct2array(load('Base_draftv2','death_year_sens'));
paper = [round(summary_HR_point(:,1)/10^6,2),...
    round(summary_HR_point(:,2)/10^6,1),...
    round(sum(inc_HR_point(11:23,:),1)/1000,1)',...
    round(summary_HR_point(:,9),0),...
    round(summary_HR_point(:,5),0),...
    round(summary_HR_point(:,8),0),...
    round(100*(summary_HR_point(:,5)-inc_HR_point(8,:))/inc_HR_point(8,:),0),...
    round(100*(summary_HR_point(:,6)-repmat(death_year_sens_point(1,6,1),length(summary_HR_point(:,6)),1))./repmat(death_year_sens_point(1,6,1),length(summary_HR_point(:,6)),1),0)];
    
% paper = [round(prctile(summary_HR(:,1,:),50,3)/10^6,2),...
%     round(prctile(summary_HR(:,2,:),50,3)/10^6,1),...
%     round(prctile(sum(inc_HR(11:23,:,:),1),50,3)/1000,1)',...
%     round(prctile(summary_HR(:,9,:),50,3),0),...
%     round(prctile(summary_HR(:,5,:),50,3),0),...
%     round(prctile(summary_HR(:,8,:),50,3),0),...
%     round(100*(prctile(summary_HR(:,5,:),50,3)-prctile(inc_HR(8,:,:),50,3))/prctile(inc_HR(8,:,:),50,3),0),...
%     round(100*(prctile(summary_HR(:,6,:),50,3)-repmat(prctile(death_year_sens(1,6,:),50,3),length(prctile(summary_HR(:,6,:),50,3)),1))./repmat(prctile(death_year_sens(1,6,:),50,3),length(prctile(summary_HR(:,6,:),50,3)),1),0)];
% 

paper_LB = [round(prctile(summary_HR(:,1,:),5,3)/10^6,2),...
    round(prctile(summary_HR(:,2,:),5,3)/10^6,1),...
    round(prctile(sum(inc_HR(11:23,:,:),1),5,3)/1000,1)',...
    round(prctile(summary_HR(:,9,:),5,3),0),...
    round(prctile(summary_HR(:,5,:),5,3),0),...
    round(prctile(summary_HR(:,8,:),5,3),0),...
    round(100*(prctile(summary_HR(:,5,:),5,3)-prctile(inc_HR(8,:,:),5,3))/prctile(inc_HR(8,:,:),5,3),0),...
    round(100*(prctile(summary_HR(:,6,:),5,3)-repmat(prctile(death_year_sens(1,6,:),5,3),length(prctile(summary_HR(:,6,:),5,3)),1))./repmat(prctile(death_year_sens(1,6,:),5,3),length(prctile(summary_HR(:,6,:),5,3)),1),0)];


paper_UB = [round(prctile(summary_HR(:,1,:),95,3)/10^6,2),...
    round(prctile(summary_HR(:,2,:),95,3)/10^6,1),...
    round(prctile(sum(inc_HR(11:23,:,:),1),95,3)/1000,1)',...
    round(prctile(summary_HR(:,9,:),95,3),0),...
    round(prctile(summary_HR(:,5,:),95,3),0),...
    round(prctile(summary_HR(:,8,:),95,3),0),...
    round(100*(prctile(summary_HR(:,5,:),95,3)-prctile(inc_HR(8,:,:),95,3))/prctile(inc_HR(8,:,:),95,3),0),...
    round(100*(prctile(summary_HR(:,6,:),95,3)-repmat(prctile(death_year_sens(1,6,:),95,3),length(prctile(summary_HR(:,6,:),95,3)),1))./repmat(prctile(death_year_sens(1,6,:),95,3),length(prctile(summary_HR(:,6,:),95,3)),1),0)];



paper_text = zeros(length(paper(:,1))+1,1+length(paper(1,:)));
paper_text = num2cell(paper_text);
paper_text{1,1} = 'Scenario';
paper_text{1,2} = 'Total costs 2018-2030 (US$M)';
paper_text{1,3} = 'Total QALYs 2018-2030 (M)';
paper_text{1,4} = 'Cumulative incidence 2018-2030 (thousand)';
paper_text{1,5} = 'Cumulative liver-related deaths 2018-2030';
paper_text{1,6} = 'Incidence in 2030';
paper_text{1,7} = 'Prevalence among PWID 2030';
paper_text{1,8} = 'Incidence reduction relative to 2015';
paper_text{1,9} = 'Mortality reduction relative to 2015';



for i = 1:length(paper(:,1))
    for j = 1:length(paper(1,:))
        paper_text{i+1,j+1} = [num2str(paper(i,j)),' (',num2str(paper_LB(i,j)),', ',num2str(paper_UB(i,j)),')'];
    end
end

if isempty(find(diag_test(:,2,2,1)>90,1))==1 aa1=0; else aa1=find(diag_test(:,2,2,1)>90,1); end
if isempty(find(diag_serum(:,2,2,1)>90,1))==1 aa2=0; else aa2=find(diag_serum(:,2,2,1)>90,1); end
if isempty(find(diag_DBS(:,2,2,1)>90,1))==1 aa3=0; else aa3=find(diag_DBS(:,2,2,1)>90,1); end

paper2_0HR = [round([0, diag_test(23,2,2,1), diag_serum(23,2,2,1), diag_DBS(23,2,2,1)],0);...
     [2007, aa1+2007,aa2+2007,aa3+2007];...
     round([summary_HR(1,1,1), summary_test(2,2,1,1,1), summary_serum(2,2,1,1,1), summary_DBS(2,2,1,1,1)]/10^6,2);...% Total costs
     round([summary_HR(1,4,1), summary_test(2,2,1,4,1), summary_serum(2,2,1,4,1), summary_DBS(2,2,1,4,1)]/10^3,2);... % total treatments
     round([summary_HR(1,9,1), summary_test(2,2,1,9,1), summary_serum(2,2,1,9,1), summary_DBS(2,2,1,9,1)],0);... % total deaths
     round([summary_HR(1,2,1), summary_test(2,2,1,2,1), summary_serum(2,2,1,2,1), summary_DBS(2,2,1,2,1)]/10^3,0);... % total DALYs
     round([summary_HR(1,5,1), summary_test(2,2,1,5,1), summary_serum(2,2,1,5,1), summary_DBS(2,2,1,5,1)],0);... % incidence in 2030
     round([summary_HR(1,8,1), summary_test(2,2,1,8,1), summary_serum(2,2,1,8,1), summary_DBS(2,2,1,8,1)],0);... % prevalence in 2030
     round(100*(inc2017 - [summary_HR(1,5,1), summary_test(2,2,1,5,1), summary_serum(2,2,1,5,1), summary_DBS(2,2,1,5,1)])/inc2017,0);... % incidence reduction in 2030
     round(100*(death2017 - [summary_HR(1,6,1), summary_test(2,2,1,6,1), summary_serum(2,2,1,6,1), summary_DBS(2,2,1,6,1)])/death2017,0)]; % mortality reduction in 2030
paper2_currentHR = [round([0, diag_test(23,2,2,2), diag_serum(23,2,2,2), diag_DBS(23,2,2,2)],0);...
     [2007, aa1+2007,aa2+2007,aa3+2007];...
     round([summary_HR(2,1,1), summary_test(2,2,2,1,1), summary_serum(2,2,2,1,1), summary_DBS(2,2,2,1,1)]/10^6,2);...% Total costs
     round([summary_HR(2,4,1), summary_test(2,2,2,4,1), summary_serum(2,2,2,4,1), summary_DBS(2,2,2,4,1)]/10^3,2);... % total treatments
     round([summary_HR(2,9,1), summary_test(2,2,2,9,1), summary_serum(2,2,2,9,1), summary_DBS(2,2,2,9,1)],0);... % total deaths
     round([summary_HR(2,2,1), summary_test(2,2,2,2,1), summary_serum(2,2,2,2,1), summary_DBS(2,2,2,2,1)]/10^3,0);... % total DALYs
     round([summary_HR(2,5,1), summary_test(2,2,2,5,1), summary_serum(2,2,2,5,1), summary_DBS(2,2,2,5,1)],0);... % incidence in 2030
     round([summary_HR(2,8,1), summary_test(2,2,2,8,1), summary_serum(2,2,2,8,1), summary_DBS(2,2,2,8,1)],0);... % prevalence in 2030
     round(100*(inc2017 - [summary_HR(2,5,1), summary_test(2,2,2,5,1), summary_serum(2,2,2,5,1), summary_DBS(2,2,2,5,1)])/inc2017,0);... % incidence reduction in 2030
     round(100*(death2017 - [summary_HR(2,6,1), summary_test(2,2,2,6,1), summary_serum(2,2,2,6,1), summary_DBS(2,2,2,6,1)])/death2017,0)]; % mortality reduction in 2030
paper2_10HR = [round([0, diag_test(23,2,2,3), diag_serum(23,2,2,3), diag_DBS(23,2,2,3)],0);...
     [2007, aa1+2007,aa2+2007,aa3+2007];...
     round([summary_HR(3,1,1), summary_test(2,2,3,1,1), summary_serum(2,2,3,1,1), summary_DBS(2,2,3,1,1)]/10^6,2);...% Total costs
     round([summary_HR(3,4,1), summary_test(2,2,3,4,1), summary_serum(2,2,3,4,1), summary_DBS(2,2,3,4,1)]/10^3,2);... % total treatments
     round([summary_HR(3,9,1), summary_test(2,2,3,9,1), summary_serum(2,2,3,9,1), summary_DBS(2,2,3,9,1)],0);... % total deaths
     round([summary_HR(3,2,1), summary_test(2,2,3,2,1), summary_serum(2,2,3,2,1), summary_DBS(2,2,3,2,1)]/10^3,0);... % total DALYs
     round([summary_HR(3,5,1), summary_test(2,2,3,5,1), summary_serum(2,2,3,5,1), summary_DBS(2,2,3,5,1)],0);... % incidence in 2030
     round([summary_HR(3,8,1), summary_test(2,2,3,8,1), summary_serum(2,2,3,8,1), summary_DBS(2,2,3,8,1)],0);... % prevalence in 2030
     round(100*(inc2017 - [summary_HR(3,5,1), summary_test(2,2,3,5,1), summary_serum(2,2,3,5,1), summary_DBS(2,2,3,5,1)])/inc2017,0);... % incidence reduction in 2030
     round(100*(death2017 - [summary_HR(3,6,1), summary_test(2,2,3,6,1), summary_serum(2,2,3,6,1), summary_DBS(2,2,3,6,1)])/death2017,0)]; % mortality reduction in 2030
paper2_20HR = [round([0, diag_test(23,2,2,4), diag_serum(23,2,2,4), diag_DBS(23,2,2,4)],0);...
     [2007, aa1+2007,aa2+2007,aa3+2007];...
     round([summary_HR(4,1,1), summary_test(2,2,4,1,1), summary_serum(2,2,4,1,1), summary_DBS(2,2,4,1,1)]/10^6,2);...% Total costs
     round([summary_HR(4,4,1), summary_test(2,2,4,4,1), summary_serum(2,2,4,4,1), summary_DBS(2,2,4,4,1)]/10^3,2);... % total treatments
     round([summary_HR(4,9,1), summary_test(2,2,4,9,1), summary_serum(2,2,4,9,1), summary_DBS(2,2,4,9,1)],0);... % total deaths
     round([summary_HR(4,2,1), summary_test(2,2,4,2,1), summary_serum(2,2,4,2,1), summary_DBS(2,2,4,2,1)]/10^3,0);... % total DALYs
     round([summary_HR(4,5,1), summary_test(2,2,4,5,1), summary_serum(2,2,4,5,1), summary_DBS(2,2,4,5,1)],0);... % incidence in 2030
     round([summary_HR(4,8,1), summary_test(2,2,4,8,1), summary_serum(2,2,4,8,1), summary_DBS(2,2,4,8,1)],0);... % prevalence in 2030
     round(100*(inc2017 - [summary_HR(4,5,1), summary_test(2,2,4,5,1), summary_serum(2,2,4,5,1), summary_DBS(2,2,4,5,1)])/inc2017,0);... % incidence reduction in 2030
     round(100*(death2017 - [summary_HR(4,6,1), summary_test(2,2,4,6,1), summary_serum(2,2,4,6,1), summary_DBS(2,2,4,6,1)])/death2017,0)]; % mortality reduction in 2030
 paper2_30HR = [round([0, diag_test(23,2,2,5), diag_serum(23,2,2,5), diag_DBS(23,2,2,5)],0);...
     [2007, aa1+2007,aa2+2007,aa3+2007];...
     round([summary_HR(5,1,1), summary_test(2,2,5,1,1), summary_serum(2,2,5,1,1), summary_DBS(2,2,5,1,1)]/10^6,2);...% Total costs
     round([summary_HR(5,4,1), summary_test(2,2,5,4,1), summary_serum(2,2,5,4,1), summary_DBS(2,2,5,4,1)]/10^3,2);... % total treatments
     round([summary_HR(5,9,1), summary_test(2,2,5,9,1), summary_serum(2,2,5,9,1), summary_DBS(2,2,5,9,1)],0);... % total deaths
     round([summary_HR(5,2,1), summary_test(2,2,5,2,1), summary_serum(2,2,5,2,1), summary_DBS(2,2,5,2,1)]/10^3,0);... % total DALYs
     round([summary_HR(5,5,1), summary_test(2,2,5,5,1), summary_serum(2,2,5,5,1), summary_DBS(2,2,5,5,1)],0);... % incidence in 2030
     round([summary_HR(5,8,1), summary_test(2,2,5,8,1), summary_serum(2,2,5,8,1), summary_DBS(2,2,5,8,1)],0);... % prevalence in 2030
     round(100*(inc2017 - [summary_HR(5,5,1), summary_test(2,2,5,5,1), summary_serum(2,2,5,5,1), summary_DBS(2,2,5,5,1)])/inc2017,0);... % incidence reduction in 2030
     round(100*(death2017 - [summary_HR(5,6,1), summary_test(2,2,5,6,1), summary_serum(2,2,5,6,1), summary_DBS(2,2,5,6,1)])/death2017,0)]; % mortality reduction in 2030
 paper2_40HR = [round([0, diag_test(23,2,2,6), diag_serum(23,2,2,6), diag_DBS(23,2,2,6)],0);...
     [2007, aa1+2007,aa2+2007,aa3+2007];...
     round([summary_HR(6,1,1), summary_test(2,2,6,1,1), summary_serum(2,2,6,1,1), summary_DBS(2,2,6,1,1)]/10^6,2);...% Total costs
     round([summary_HR(6,4,1), summary_test(2,2,6,4,1), summary_serum(2,2,6,4,1), summary_DBS(2,2,6,4,1)]/10^3,2);... % total treatments
     round([summary_HR(6,9,1), summary_test(2,2,6,9,1), summary_serum(2,2,6,9,1), summary_DBS(2,2,6,9,1)],0);... % total deaths
     round([summary_HR(6,2,1), summary_test(2,2,6,2,1), summary_serum(2,2,6,2,1), summary_DBS(2,2,6,2,1)]/10^3,0);... % total DALYs
     round([summary_HR(6,5,1), summary_test(2,2,6,5,1), summary_serum(2,2,6,5,1), summary_DBS(2,2,6,5,1)],0);... % incidence in 2030
     round([summary_HR(6,8,1), summary_test(2,2,6,8,1), summary_serum(2,2,6,8,1), summary_DBS(2,2,6,8,1)],0);... % prevalence in 2030
     round(100*(inc2017 - [summary_HR(6,5,1), summary_test(2,2,6,5,1), summary_serum(2,2,6,5,1), summary_DBS(2,2,6,5,1)])/inc2017,0);... % incidence reduction in 2030
     round(100*(death2017 - [summary_HR(6,6,1), summary_test(2,2,6,6,1), summary_serum(2,2,6,6,1), summary_DBS(2,2,6,6,1)])/death2017,0)]; % mortality reduction in 2030
 paper2_50HR = [round([0, diag_test(23,2,2,7), diag_serum(23,2,2,7), diag_DBS(23,2,2,7)],0);...
     [2007, aa1+2007,aa2+2007,aa3+2007];...
     round([summary_HR(7,1,1), summary_test(2,2,7,1,1), summary_serum(2,2,7,1,1), summary_DBS(2,2,7,1,1)]/10^6,2);...% Total costs
     round([summary_HR(7,4,1), summary_test(2,2,7,4,1), summary_serum(2,2,7,4,1), summary_DBS(2,2,7,4,1)]/10^3,2);... % total treatments
     round([summary_HR(7,9,1), summary_test(2,2,7,9,1), summary_serum(2,2,7,9,1), summary_DBS(2,2,7,9,1)],0);... % total deaths
     round([summary_HR(7,2,1), summary_test(2,2,7,2,1), summary_serum(2,2,7,2,1), summary_DBS(2,2,7,2,1)]/10^3,0);... % total DALYs
     round([summary_HR(7,5,1), summary_test(2,2,7,5,1), summary_serum(2,2,7,5,1), summary_DBS(2,2,7,5,1)],0);... % incidence in 2030
     round([summary_HR(7,8,1), summary_test(2,2,7,8,1), summary_serum(2,2,7,8,1), summary_DBS(2,2,7,8,1)],0);... % prevalence in 2030
     round(100*(inc2017 - [summary_HR(7,5,1), summary_test(2,2,7,5,1), summary_serum(2,2,7,5,1), summary_DBS(2,2,7,5,1)])/inc2017,0);... % incidence reduction in 2030
     round(100*(death2017 - [summary_HR(7,6,1), summary_test(2,2,7,6,1), summary_serum(2,2,7,6,1), summary_DBS(2,2,7,6,1)])/death2017,0)]; % mortality reduction in 2030

paper2 = [paper2_0HR; zeros(1,4); paper2_currentHR; zeros(1,4); paper2_10HR; zeros(1,4); paper2_20HR; zeros(1,4); paper2_30HR; zeros(1,4); paper2_40HR; zeros(1,4); paper2_50HR];
 
paper2_text = zeros(length(paper2(:,1))+7,1+length(paper2(1,:)));
paper2_text = num2cell(paper2_text);
paper2_text{1,1} = 'Scenario';
paper2_text{1,2} = 'No testting/treatment';
paper2_text{1,3} = 'Ab+RNA';
paper2_text{1,4} = 'Serum-based HCVcAg';
paper2_text{1,5} = 'Dry Blood Spot HCVcAg';

paper2_text{2,1} = 'Percentage diagnosed in 2030';
paper2_text{3,1} = 'Year 80% diagnosed reached';
paper2_text{4,1} = 'Total costs 2018-2030 (million US$)';
paper2_text{5,1} = 'Total treatments 2018-2030 (thousand)';
paper2_text{6,1} = 'Total HCV-related deaths 2018-2030';
paper2_text{7,1} = 'Total DALYs 2018-2030 (million)';
paper2_text{8,1} = 'Incidence in 2030';
paper2_text{9,1} = 'Prevalence in 2030';
paper2_text{10,1} = 'Incidence reduction in 2030';
paper2_text{11,1} = 'Mortality reduction in 2030';

paper2_text{13,1} = 'Percentage diagnosed in 2030';
paper2_text{14,1} = 'Year 80% diagnosed reached';
paper2_text{15,1} = 'Total costs 2018-2030 (million US$)';
paper2_text{16,1} = 'Total treatments 2018-2030 (thousand)';
paper2_text{17,1} = 'Total HCV-related deaths 2018-2030';
paper2_text{18,1} = 'Total DALYs 2018-2030';
paper2_text{19,1} = 'Incidence in 2030';
paper2_text{20,1} = 'Prevalence in 2030';
paper2_text{21,1} = 'Incidence reduction in 2030';
paper2_text{22,1} = 'Mortality reduction in 2030';

paper2_text{24,1} = 'Percentage diagnosed in 2030';
paper2_text{25,1} = 'Year 80% diagnosed reached';
paper2_text{26,1} = 'Total costs 2018-2030 (million US$)';
paper2_text{27,1} = 'Total treatments 2018-2030 (thousand)';
paper2_text{28,1} = 'Total HCV-related deaths 2018-2030';
paper2_text{29,1} = 'Total DALYs 2018-2030';
paper2_text{30,1} = 'Incidence in 2030';
paper2_text{31,1} = 'Prevalence in 2030';
paper2_text{32,1} = 'Incidence reduction in 2030';
paper2_text{33,1} = 'Mortality reduction in 2030';

paper2_text{35,1} = 'Percentage diagnosed in 2030';
paper2_text{36,1} = 'Year 80% diagnosed reached';
paper2_text{37,1} = 'Total costs 2018-2030 (million US$)';
paper2_text{38,1} = 'Total treatments 2018-2030 (thousand)';
paper2_text{39,1} = 'Total HCV-related deaths 2018-2030';
paper2_text{40,1} = 'Total DALYs 2018-2030';
paper2_text{41,1} = 'Incidence in 2030';
paper2_text{42,1} = 'Prevalence in 2030';
paper2_text{43,1} = 'Incidence reduction in 2030';
paper2_text{44,1} = 'Mortality reduction in 2030';

paper2_text{46,1} = 'Percentage diagnosed in 2030';
paper2_text{47,1} = 'Year 80% diagnosed reached';
paper2_text{48,1} = 'Total costs 2018-2030 (million US$)';
paper2_text{49,1} = 'Total treatments 2018-2030 (thousand)';
paper2_text{50,1} = 'Total HCV-related deaths 2018-2030';
paper2_text{51,1} = 'Total DALYs 2018-2030';
paper2_text{52,1} = 'Incidence in 2030';
paper2_text{53,1} = 'Prevalence in 2030';
paper2_text{54,1} = 'Incidence reduction in 2030';
paper2_text{55,1} = 'Mortality reduction in 2030';

paper2_text{57,1} = 'Percentage diagnosed in 2030';
paper2_text{58,1} = 'Year 80% diagnosed reached';
paper2_text{59,1} = 'Total costs 2018-2030 (million US$)';
paper2_text{60,1} = 'Total treatments 2018-2030 (thousand)';
paper2_text{61,1} = 'Total HCV-related deaths 2018-2030';
paper2_text{62,1} = 'Total DALYs 2018-2030';
paper2_text{63,1} = 'Incidence in 2030';
paper2_text{64,1} = 'Prevalence in 2030';
paper2_text{65,1} = 'Incidence reduction in 2030';
paper2_text{66,1} = 'Mortality reduction in 2030';

paper2_text{68,1} = 'Percentage diagnosed in 2030';
paper2_text{69,1} = 'Year 80% diagnosed reached';
paper2_text{70,1} = 'Total costs 2018-2030 (million US$)';
paper2_text{71,1} = 'Total treatments 2018-2030 (thousand)';
paper2_text{72,1} = 'Total HCV-related deaths 2018-2030';
paper2_text{73,1} = 'Total DALYs 2018-2030';
paper2_text{74,1} = 'Incidence in 2030';
paper2_text{75,1} = 'Prevalence in 2030';
paper2_text{76,1} = 'Incidence reduction in 2030';
paper2_text{77,1} = 'Mortality reduction in 2030';

for i = 1:length(paper2(:,1))
    for j = 1:length(paper2(1,:))
        paper2_text{i+1,j+1} = [num2str(paper2(i,j))];%,' (',num2str(paper_LB(i,j)),', ',num2str(paper_UB(i,j)),')'];
    end
end


% paper_sens = [round(summary_DBS_rr(1,:,3,1)'/10^6,2),...
%     round(summary_DBS_rr(1,:,3,2)/10^6,1)',...
%     round(reshape(sum(inc_DBS_rr(11:23,1,:,3),1)/1000,length(range_diagnosed_risk_reduction),1),1),...
%     round(summary_DBS_rr(1,:,3,9),0)',...
%     round(summary_DBS_rr(1,:,3,5),0)',...
%     round(summary_DBS_rr(1,:,3,8),0)',...
%     round(100*(summary_DBS_rr(1,:,3,5)'-reshape(inc_DBS_rr(8,1,:,3),length(range_diagnosed_risk_reduction),1))...
%     ./reshape(inc_DBS_rr(8,1,:,3),length(range_diagnosed_risk_reduction),1),0),...
%     round(100*(summary_DBS_rr(1,:,3,6)'-repmat(death_year_sens(1,6,1),length(range_diagnosed_risk_reduction),1))./repmat(death_year_sens(1,6,1),length(range_diagnosed_risk_reduction),1),0)];
%  
% paper3_text = zeros(length(paper_sens(:,1))+1,1+length(paper_sens(1,:)));
% paper3_text = num2cell(paper3_text);
% paper3_text{1,1} = 'Scenario';
% paper3_text{1,2} = 'Total costs 2018-2030 (US$M)';
% paper3_text{1,3} = 'Total QALYs 2018-2030 (M)';
% paper3_text{1,4} = 'Cumulative incidence 2018-2030 (thousand)';
% paper3_text{1,5} = 'Cumulative liver-related deaths 2018-2030';
% paper3_text{1,6} = 'Incidence in 2030';
% paper3_text{1,7} = 'Prevalence among PWID 2030';
% paper3_text{1,8} = 'Incidence reduction relative to 2015';
% paper3_text{1,9} = 'Mortality reduction relative to 2015';
% 
% 
% for i = 1:length(paper_sens(:,1))
%     for j = 1:length(paper_sens(1,:))
%         paper3_text{i+1,j+1} = [num2str(paper_sens(i,j)),' (',num2str(paper_LB(i,j)),', ',num2str(paper_UB(i,j)),')'];
%     end
% end

