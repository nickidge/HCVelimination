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

%Define the output file using variables based on the current
%working directory or relative paths so that this can work on
%either Nick or Rachel's computer - this will only work if you are in a
%subdirectory of ..\Users\User not in the shared drives
%filename = 'C:\Users\Nick\Desktop\Matlab Sims\Testing\foo';
user=extractBetween(pwd,"Users\","\");
drive=extractBefore(pwd,":");
filename=strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\foo5");
%relative path method is commented out below because it's less general than
%the above method.
%filename='../../../Desktop/Matlab Sims/Testing/foo';

loaddata
dt = 1/4; % six-monthly time steps for burn-in / calibration perio 1950-2015
sens=1; %Number of runs in sensitivity analysis, sens=1 turns off feature
summary=zeros(6,12,sens);
followup = 1;
    

%% Sensitivity of parameters
for s=1:sens
    s
    scenario = 'empty';
    if sens > 1
        perturb_parameters
    end
    
    alpha=alpha_old; % calibrate in pre-DAA era
    target_late=1; % target treatments to people with late-liver disease

    %moved these up because need to use these for calibration if
    %calibration is running.
    infect_base=infect;
    progression_base = progression;
    
    
    [output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases,output_ost,output_nsp, output_diagnoses] = ...
        calibrate_optim_par(200, 30);
    save('calibration_test')
    %load('C:\Users\Nick\Desktop\Matlab Sims\Tanzania\calibration_draftv1')
    %filename is stored in calibration_data so have to add here so that we
    %can have multiple users using these files
    %filename=strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\foo");
    
    infect_base=infect;
    progression_base = progression;
    dt = 1/4;
    %Run model in prior to 2017
    [TT1,y1]=DE_track_age(Tin,y0,t0,treat);
    y1_end=reshape(y1(end,:,:,:,:,:,:,:), num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region,33);
    death2017 = (sum(sum(sum(sum(sum(sum(y1(end,:,:,:,:,:,:,22))))))) - sum(sum(sum(sum(sum(sum(y1(find(TT1>=66,1),:,:,:,:,:,:,22)))))))) / ...
        (TT1(end)-TT1(find(TT1>=66,1)));
    inc2017 = sum(sum(sum(sum(sum(sum(sum(y1(find(TT1>=66,1):end-1,:,:,:,:,:,:,27))))))));%/(TT1(find(TT1>=66,1)-1) - TT1(find(TT1>=65,1)));
    %y1_end(:,:,:,:,:,:,21:(27+6))=0;
    y1_end(:,6,:,:,:,:,1:20) = y1_end(:,6,:,:,:,:,1:20) + y1_end(:,8,:,:,:,:,1:20) + y1_end(:,10,:,:,:,:,1:20); 
    y1_end(:,8,:,:,:,:,1:20) = 0; y1_end(:,10,:,:,:,:,1:20) = 0; % moving failed to be re-eligible with DAAs
    
    sum(sum(sum(sum(sum(y1(find(TT1>=Tin-(2016-prev0(1,1)),1),1,:,:,:,:,1,6:20))))))./sum(sum(sum(sum(sum(y1(find(TT1>=Tin-(2016-prev0(1,1)),1),1,:,:,:,:,1,1:20))))))
    
    %% Baseline: Current standard of care with no scaled up treatment
    scenario = 'base'; %Current level of community care
    alpha = alpha_DAA;
    prop_test = 1;
     
    [TT2,y2]=DE_track_age(Run,y1_end,TT1,treat);
    [ycomb_noage, summary(1,:,s), tr, tr_] = gather_outputs(y1,y2,TT2);
    
 
    %%  Scenario 1: Scaled harm reduction
    scenario = 'current'
    harm_reduction_range = [0,0.06,0.1:0.1:0.5];
    summary_HR = zeros(length(harm_reduction_range),length(summary(1,:,1)),sens);
    for h = 1:length(harm_reduction_range)
        nsp_coverage = harm_reduction_range(h);
        ost_coverage = harm_reduction_range(h);
        [TT2_treat,y2_treat]=DE_track_age(Run,y1_end,TT1,treat);
        [ycomb2_noage, summary(2,:,s), tr2, tr2_] = gather_outputs(y1,y2_treat,TT2_treat);
        %         for i = Tin-10:Tin
        %             prev_HR(i,h) = 100*sum(sum(sum(sum(sum(sum(y1(find(TT1>=i,1),1,:,:,:,:,:,[6,12:20])))))))./...
        %                 sum(sum(sum(sum(sum(sum(y1(find(TT1>=i,1),1,:,:,:,:,:,1:20)))))));
        %             inc_HR(i,h) = sum(sum(sum(sum(sum(sum(sum(y1(find(TT1>=i-1,1):find(TT1>=i,1)-1,:,:,:,:,:,:,27))))))));
        %         end
        for i = (Tin-10):(Tin+Run)
            inc_HR(i-Tin+11,h) = (sum(sum(sum(ycomb2_noage(find(TT2_treat>=i-1,1):find(TT2_treat>=i,1)-1,:,:,27)))));
            prev_HR(i-Tin+11,h) = 100*sum(sum(ycomb2_noage(find(TT2_treat>=i,1),1,:,[6,12:20])))./...
                sum(sum(ycomb2_noage(find(TT2_treat>=i,1),1,:,1:20)));
        end
        summary_HR(h,:,s) = summary(2,:,s);
    end
    nsp_coverage = 0.06;
    ost_coverage = 0.06;
    
    %%  Scenario 2: Ab + RNA (standard testing) to reach 90% diagnosed
    scenario = 'current'; dt = 1/4;
    target_late=0; % Target PWID
    alpha = alpha_DAA;
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    range_test = [0.5, 1, 2, 4]*2; % divided by 2 to assume that infection happens midway between tests
    range_followup = [0.74];
    prop_test_range = [0.6,0.8,0.9];
    for i = 1:length(range_test)
        for j = 1:length(range_followup)
            for k =1:length(prop_test_range)
                prop_test = prop_test_range(k);
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
                if range_followup(j) > 0
                    followup = range_followup(j);%/0.46;
                    progression(2,2,2,1) = 1/0.25; progression(1,2,2,1) = 1 / 0.25; progression(3,2,2,1) = 1 / 0.25;
                else
                    progression(1,2,2,1) = progression_base(1,2,2,1);
                    progression(2,2,2,1) = progression_base(2,2,2,1);
                    progression(3,2,2,1) = progression_base(3,2,2,1);
                    followup = 1;
                end
                y1_end_sim = y1_end;
                y1_end_sim(:,:,:,:,1,:,:) = (1-prop_test)* y1_end(:,:,:,:,2,:,:);
                y1_end_sim(:,:,:,:,2,:,:) = prop_test * y1_end(:,:,:,:,2,:,:);
                y1_end_sim(:,1,:,:,1,:,1:20) = y1_end_sim(:,1,:,:,1,:,1:20) + sum(y1_end_sim(:,2:end,:,:,1,:,1:20),2);
                y1_end_sim(:,2:end,:,:,1,:,1:20) = 0;
                [TT2_treat5,y2_treat5]=DE_track_age(Run,y1_end_sim,TT1,treat);
                [ycomb5_noage, summary(5,:,s), tr5, tr5_] = gather_outputs(y1,y2_treat5,TT2_treat5);
                
                for t = (Tin-10):(Tin+Run)
                    inc_test(t-Tin+11,i,j,k) = (sum(sum(sum(ycomb5_noage(find(TT2_treat5>=t-1,1):find(TT2_treat5>=t,1)-1,:,:,27)))));
                    prev_test(t-Tin+11,i,j,k) = 100*sum(sum(ycomb5_noage(find(TT2_treat5>=t,1),1,:,[6,12:20])))./...
                        sum(sum(ycomb5_noage(find(TT2_treat5>=t,1),1,:,1:20)));
                    diag_test(t-Tin+11,i,j,k) = 100*sum(sum(sum(ycomb5_noage(find(TT2_treat5>=t,1),:,3:end,[6,12:20]))))./...
                        sum(sum(sum(ycomb5_noage(find(TT2_treat5>=t,1),:,:,[6,12:20]))));
                end
                summary_test(i,j,k,:,s) = summary(5,:,s);
            end
        end
    end
    
    %%  Scenario 3:  Serum-based HCVcAg testing
    scenario = 'DBS_HCVcAg';
    target_late=0; % Target PWID
    alpha = alpha_DAA;
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    progression(1,2,2,1) = 50; progression(2,2,2,1) = 50; progression(3,2,2,1) = 50; % perfect "RNA" followup
    range_test = [0.5, 1, 2, 4]*2; % divided by 2 to assume that infection happens midway between tests
    range_followup = [0.76]; % 76% sensitivity
    prop_test_range = [0.6,0.8,0.9];
    for i = 1:length(range_test)
        for j = 1:length(range_followup)
            for k =1:length(prop_test_range)
                prop_test = prop_test_range(k);
                followup = range_followup(j);
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
                [TT2_treat3,y2_treat3]=DE_track_age(Run,y1_end_sim,TT1,treat);
                [ycomb3_noage, summary(3,:,s), tr3, tr3_] = gather_outputs(y1,y2_treat3,TT2_treat3);
                
                for t = (Tin-10):(Tin+Run)
                    inc_DBS(t-Tin+11,i,j,k) = (sum(sum(sum(ycomb3_noage(find(TT2_treat3>=t-1,1):find(TT2_treat3>=t,1)-1,:,:,27)))));
                    prev_DBS(t-Tin+11,i,j,k) = 100*sum(sum(ycomb3_noage(find(TT2_treat3>=t,1),1,:,[6,12:20])))./...
                        sum(sum(ycomb3_noage(find(TT2_treat3>=t,1),1,:,1:20)));
                    diag_DBS(t-Tin+11,i,j,k) = 100*sum(sum(sum(ycomb3_noage(find(TT2_treat3>=t,1),:,2:end,[6,12:20]))))./...
                        sum(sum(sum(ycomb3_noage(find(TT2_treat3>=t,1),:,:,[6,12:20]))));
                end
                summary_DBS(i,j,k,:,s) = summary(3,:,s);
            end
        end
    end

    
    %%  Scenario 4: DBS HCVcAg testing
    scenario = 'serum_HCVcAg';
    target_late=0; % Target PWID
    alpha = alpha_DAA;
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    progression(1,2,2,1) = 50; progression(2,2,2,1) = 50; progression(3,2,2,1) = 50; % perfect "RNA" followup
    range_test = [0.5, 1, 2, 4]*2; % divided by 2 to assume that infection happens midway between tests
    range_followup = [1]; % 100% sensitivity
    prop_test_range = [0.6,0.8,0.9];
    for i = 1:length(range_test)
        for j = 1:length(range_followup)
            for k =1:length(prop_test_range)
                prop_test = prop_test_range(k);
                followup = range_followup(j);
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
                [TT2_treat6,y2_treat6]=DE_track_age(Run,y1_end_sim,TT1,treat);
                [ycomb6_noage, summary(6,:,s), tr6, tr6_] = gather_outputs(y1,y2_treat6,TT2_treat6);
                
                for t = (Tin-10):(Tin+Run)
                    inc_serum(t-Tin+11,i,j,k) = (sum(sum(sum(ycomb6_noage(find(TT2_treat6>=t-1,1):find(TT2_treat6>=t,1)-1,:,:,27)))));
                    prev_serum(t-Tin+11,i,j,k) = 100*sum(sum(ycomb6_noage(find(TT2_treat6>=t,1),1,:,[6,12:20])))./...
                        sum(sum(ycomb6_noage(find(TT2_treat6>=t,1),1,:,1:20)));
                    diag_serum(t-Tin+11,i,j,k) = 100*sum(sum(sum(ycomb6_noage(find(TT2_treat6>=t,1),:,2:end,[6,12:20]))))./...
                        sum(sum(sum(ycomb6_noage(find(TT2_treat6>=t,1),:,:,[6,12:20]))));
                end
                summary_serum(i,j,k,:,s) = summary(6,:,s);
            end
        end
    end
    
    
    %%  Scenario 5: Diagnosed risk reduction
    scenario = 'DBC_HCVcAg';
    target_late=0; % Target PWID
    alpha = alpha_DAA;
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    progression(1,2,2,1) = 50; progression(2,2,2,1) = 50; progression(3,2,2,1) = 50; % perfect "RNA" followup
    range_test = [0.5, 1, 2, 4]*2; % divided by 2 to assume that infection happens midway between tests
    range_diagnosed_risk_reduction = [0:0.1:0.5]; % risk reduction when diagnosed
    prop_test_range = [0.6,0.8,0.9];
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


inc_year=mean(inc_year_sens(:,:,:)./repmat(inc_year_sens(1,1,:),6,length(years)-1),3);
%inc_year(5,:)=mean(inc_year_sens(5,:,treat_scaleup(1,:)>=0)./repmat(inc_year_sens(1,1,treat_scaleup(1,:)>=0),1,length(years)-1),3);

%inc_year = mean(inc_year_sens(:,:,:),3);
inc_year_std=prctile((inc_year_sens(:,:,:)./repmat(inc_year_sens(1,1,:),6,length(years)-1)),[5,95],3);
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
Charts
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


paper = [round(summary_HR(:,1)/10^6,2),...
    round(summary_HR(:,2)/10^6,1),...
    round(sum(inc_HR(11:23,:),1)/1000,1)',...
    round(summary_HR(:,9),0),...
    round(summary_HR(:,5),0),...
    round(summary_HR(:,8),0),...
    round(100*(summary_HR(:,5)-inc_HR(8,:))/inc_HR(8,:),0),...
    round(100*(summary_HR(:,6)-repmat(death_year_sens(1,6,1),length(summary_HR(:,6)),1))./repmat(death_year_sens(1,6,1),length(summary_HR(:,6)),1),0)];
    
    
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
        paper_text{i+1,j+1} = [num2str(paper(i,j))];%,' (',num2str(paper_LB(i,j)),', ',num2str(paper_UB(i,j)),')'];
    end
end

if isempty(find(diag_test(:,1,1,2)>90,1))==1 aa=0; else aa=find(diag_test(:,1,1,2)>90,1); end

paper2 = [round([diag_test(23,1,1,2); diag_serum(23,1,1,2); diag_DBS(23,1,1,3)],0),...
    [aa+2007;find(diag_serum(:,1,1,2)>90,1)+2007;find(diag_DBS(:,1,1,3)>90,1)+2007],...
    round([summary_test(1,1,2,1,1); summary_serum(1,1,2,1,1); summary_DBS(1,1,3,1,1)]/10^6,2)];
paper2_text = zeros(length(paper2(:,1))+1,1+length(paper2(1,:)));
paper2_text = num2cell(paper2_text);
paper2_text{1,1} = 'Scenario';
paper2_text{1,2} = 'Percentage diagnosed in 2030';
paper2_text{1,3} = 'Year 80% diagnosis target reached';
paper2_text{1,4} = 'Total cost';

for i = 1:length(paper2(:,1))
    for j = 1:length(paper2(1,:))
        paper2_text{i+1,j+1} = [num2str(paper2(i,j))];%,' (',num2str(paper_LB(i,j)),', ',num2str(paper_UB(i,j)),')'];
    end
end


paper_sens = [round(summary_DBS_rr(1,:,3,1)'/10^6,2),...
    round(summary_DBS_rr(1,:,3,2)/10^6,1)',...
    round(reshape(sum(inc_DBS_rr(11:23,1,:,3),1)/1000,length(range_diagnosed_risk_reduction),1),1),...
    round(summary_DBS_rr(1,:,3,9),0)',...
    round(summary_DBS_rr(1,:,3,5),0)',...
    round(summary_DBS_rr(1,:,3,8),0)',...
    round(100*(summary_DBS_rr(1,:,3,5)'-reshape(inc_DBS_rr(8,1,:,3),length(range_diagnosed_risk_reduction),1))...
    ./reshape(inc_DBS_rr(8,1,:,3),length(range_diagnosed_risk_reduction),1),0),...
    round(100*(summary_DBS_rr(1,:,3,6)'-repmat(death_year_sens(1,6,1),length(range_diagnosed_risk_reduction),1))./repmat(death_year_sens(1,6,1),length(range_diagnosed_risk_reduction),1),0)];
 
paper3_text = zeros(length(paper_sens(:,1))+1,1+length(paper_sens(1,:)));
paper3_text = num2cell(paper3_text);
paper3_text{1,1} = 'Scenario';
paper3_text{1,2} = 'Total costs 2018-2030 (US$M)';
paper3_text{1,3} = 'Total QALYs 2018-2030 (M)';
paper3_text{1,4} = 'Cumulative incidence 2018-2030 (thousand)';
paper3_text{1,5} = 'Cumulative liver-related deaths 2018-2030';
paper3_text{1,6} = 'Incidence in 2030';
paper3_text{1,7} = 'Prevalence among PWID 2030';
paper3_text{1,8} = 'Incidence reduction relative to 2015';
paper3_text{1,9} = 'Mortality reduction relative to 2015';


for i = 1:length(paper_sens(:,1))
    for j = 1:length(paper_sens(1,:))
        paper3_text{i+1,j+1} = [num2str(paper_sens(i,j))];%,' (',num2str(paper_LB(i,j)),', ',num2str(paper_UB(i,j)),')'];
    end
end

