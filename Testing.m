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
    harm_reduction_coverage nsp_coverage ost_coverage

%Define the output file using variables based on the current
%working directory or relative paths so that this can work on
%either Nick or Rachel's computer - this will only work if you are in a
%subdirectory of ..\Users\User not in the shared drives
%filename = 'C:\Users\Nick\Desktop\Matlab Sims\Testing\foo';
user=extractBetween(pwd,"Users\","\");
drive=extractBefore(pwd,":");
filename=strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\foo2");
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
        calibrate_optim_par(200, 20);
    save('C:\Users\Nick\Desktop\Matlab Sims\Tanzania\calibration2')
    load('C:\Users\Nick\Desktop\Matlab Sims\Tanzania\calibration2')
    %filename is stored in calibration_data so have to add here so that we
    %can have multiple users using these files
    %filename=strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\foo");
    
    infect_base=infect;
    progression_base = progression;

    %Run model in prior to 2017
    [TT1,y1]=DE_track_age(Tin,y0,t0,treat);
    y1_end=reshape(y1(end,:,:,:,:,:,:,:), num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region,33);
    death2017 = (sum(sum(sum(sum(sum(sum(y1(end,:,:,:,:,:,:,22))))))) - sum(sum(sum(sum(sum(sum(y1(find(TT1>=66,1),:,:,:,:,:,:,22)))))))) / ...
        (TT1(end)-TT1(find(TT1>=66,1)));
    inc2017 = sum(sum(sum(sum(sum(sum(sum(y1(find(TT1>=66,1):end-1,:,:,:,:,:,:,27))))))));%/(TT1(find(TT1>=66,1)-1) - TT1(find(TT1>=65,1)));
    y1_end(:,:,:,:,:,:,21:(27+6))=0;
    y1_end(:,6,:,:,:,:,1:20) = y1_end(:,6,:,:,:,:,1:20) + y1_end(:,8,:,:,:,:,1:20) + y1_end(:,10,:,:,:,:,1:20); 
    y1_end(:,8,:,:,:,:,1:20) = 0; y1_end(:,10,:,:,:,:,1:20) = 0; % moving failed to be re-eligible with DAAs
    
    
    %% Baseline: Current standard of care with no scaled up treatment
    scenario = 'base'; %Current level of community care
    alpha = alpha_DAA;
    dt = 1/4; 
    [TT2,y2]=DE_track_age(Run,y1_end,TT1,treat);
    [ycomb_noage, summary(1,:,s), tr, tr_] = gather_outputs(y1,y2,TT2);
    
 
    %%  Scenario 1: Scaled harm reduction
    scenario = 'current'
    harm_reduction_range = [0,0.06,0.1:0.1:0.5];
    summary_HR = zeros(length(harm_reduction_range),length(summary(1,:,1)),sens);
    for h = 1:length(harm_reduction_range)
        nsp_coverage = harm_reduction_range(h)
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
                sum(sum(ycomb_noage(find(TT2_treat>=i,1),1,:,1:20)));
        end
        summary_HR(h,:,s) = summary(2,:,s);
    end

    %%  Scenario 2: Calculate required treatments for incidence target
    scenario = 'current'
    target_late=0; % Target PWID
    alpha = alpha_DAA;
    %apply_cascade_rates(scenario, 0);
    for i =1:num_pops
        for j =1:6
            for k =2:num_engagement
                progression(i,j,k,1) = 50;
            end
        end
    end

    elim_years = [80];
    if strcmp(scenario,'current') == 1
        for j = 1:length(elim_years)
            treat_range = 2000:10000:55000; % Estimate treatments required
            for i = 1:length(treat_range)
                i
                treat(1) = treat_range(i); treat(2) = 0;
                [TT2_treat5,y2_treat5]=DE_track_age(Run,y1_end,TT1,treat);
                [ycomb5_noage, summary(5,:,s), tr5, tr5_] = gather_outputs(y1,y2_treat5,TT2_treat5);
                
                incred(i,j) = (inc2017 - sum(sum(sum(ycomb5_noage(find(TT2_treat5>=elim_years(j),1):find(TT2_treat5>=elim_years(j)+1,1)-1,:,:,27)))))/ inc2017 ;%...
                    %(TT2_treat5(find(TT2_treat5>=elim_years(j)+1,1)-1)-TT2_treat5(find(TT2_treat5>=elim_years(j),1)))) / inc2015;
                if incred(i,j) > target_inc && i > 1
                    break;
                end
            end
            if incred(i,j) < target_inc
                treat_scaleup(j,s) = -1;
            else
                treat_scaleup(j,s) = interp1(incred(1:i,j), treat_range(1:i), target_inc);
                treat(1) = treat_scaleup(j,s);
            end
        end
    end
    treat(1) = treat_scaleup(1,s); treat(2) = 0;
    [TT2_treat5,y2_treat5]=DE_track_age(Run,y1_end,TT1,treat);
    [ycomb5_noage, summary(5,:,s), tr5, tr5_] = gather_outputs(y1,y2_treat5,TT2_treat5);
    
    
    %%  Scenario 3:  Eliminate by 2020
    scenario = 'current'
    target_late=0; % Target PWID
    alpha = alpha_DAA;
    for i =1:num_pops
        for j =1:6
            for k =1:num_engagement
                progression(i,j,k,1) = 50;
            end
        end
    end

    %treat(1) = treat_scaleup(2,s); treat(2) = 0;
    
    [TT2_treat6,y2_treat6]=DE_track_age(Run,y1_end,TT1,treat);
    [ycomb6_noage, summary(6,:,s), tr6, tr6_] = gather_outputs(y1,y2_treat6,TT2_treat6);
    
    
    %%  Scenario 4: antibody testing rates
    scenario = 'current';
    progression = progression_base;
    progression(1,5,2,1) = 1/(30/365); progression(2,5,2,1) = 1/(30/365); progression(3,5,2,1) = 1/(30/365);
    progression(1,6,2,1) = 1/(30/365); progression(2,6,2,1) = 1/(30/365); progression(3,6,2,1) = 1/(30/365);
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    target_late = -0.5;
    range_test = [-1, 1/2, 2, 1];
    range_followup = [-1,0:0.2:1];
    treat = treat_projected; 
    for i = 1:length(range_test)
        for j = 1:length(range_followup)
            if range_test(i) > 0
                temp_rate = (progression_base(1,1,2,1));% - 0.9*0.5) / (0.1);
                progression(1,1,2,1) = range_test(i);% + 0.1*temp_rate;
                progression(2,1,2,1) = range_test(i);% + 0.1*temp_rate;
                %progression(3,1,2,1) = 0.9*range_test(i) + 0.1*temp_rate;
            end
            if range_followup(j) > 0
                temp_rate1 = progression_base(1,2,2,1);
                temp_rate2 = progression_base(2,2,2,1);
                progression(1,2,2,1) = range_followup(j)*(1/0.25) + (1-range_followup(j))*temp_rate1;
                progression(2,2,2,1) = range_followup(j)*(1/0.25) + (1-range_followup(j))*temp_rate2;
            else
                progression(1,2,2,1) = progression_base(1,2,2,1);
                progression(2,2,2,1) = progression_base(2,2,2,1);
            end
            
            [TT2_treat3,y2_treat3]=DE_track_age(Run,y1_end,TT1,treat);
            [ycomb3_noage, summary(3,:,s), tr3, tr3_] = gather_outputs(y1,y2_treat3,TT2_treat3);
            
            for year_elim = 1:Run-2
                incred_test(i,year_elim) = (inc2017 - sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin+year_elim,1):find(TT2_treat3>=Tin+year_elim+1,1)-1,:,:,27)))))/...
                    (TT2_treat3(find(TT2_treat3>=Tin+year_elim+1,1)-1)-TT2_treat3(find(TT2_treat3>=Tin+year_elim,1))) / inc2017;
            end
            if incred_test(i,end) >= 0.8
                year_elim_test(i,j) = 1950+Tin+find(incred_test(i,:)>=0.8,1);
            else
                year_elim_test(i,j) = -1;
            end
        end
    end
    
    
    %%  Scenario 5: RNA testing
    scenario = 'annualRNA';
    progression = progression_base;
    progression(1,5,2,1) = 1/(30/365); progression(2,5,2,1) = 1/(30/365); progression(3,5,2,1) = 1/(30/365);
    progression(1,6,2,1) = 1/(30/365); progression(2,6,2,1) = 1/(30/365); progression(3,6,2,1) = 1/(30/365);
    progression(1,3,2,1) = 50; progression(2,3,2,1) = 50; progression(3,3,2,1) = 50; % remove genotype
    target_late = -0.5; 
    treat = treat_projected;
    [TT2_treat4,y2_treat4]=DE_track_age(Run,y1_end,TT1,treat);
    [ycomb4_noage, summary(4,:,s), tr4, tr4_] = gather_outputs(y1,y2_treat4,TT2_treat4);

    
   
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
    death_year_sens(:,6,s) = [death2015;death2015;death2015;death2015;death2015;death2015];
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

paper = round([total_treat_summary([1,2,5,6,4],:,1),...
    total_treat_2030_summary([1,2,5,6,4],4,1),...
    [treat_projected(3); sum(treat_projected(1:2));treat_scaleup(1,1);treat_scaleup(1,1);sum(treat_projected(1:2))],...
    [margins(1:2,1,1); margins(5:6,1,1); margins(4,1,1)],...
    [margins(1:2,2,1); margins(5:6,2,1); margins(4,2,1)],...
    [margins(1:2,[4,5,6,8,9],1); margins(5:6,[4,5,6,8,9],1); margins(4,[4,5,6,8,9],1)],...
    100-100*[inc_year(1:2,10); inc_year(5:6,10); inc_year(4,10)],... %inc reduction in 2018
    100-100*[inc_year(1:2,11); inc_year(5:6,11); inc_year(4,11)],...%inc reduction in 2020
    100-100*[inc_year(1:2,16); inc_year(5:6,16); inc_year(4,16)],...%inc reduction in 2025
    100-100*[inc_year(1:2,21); inc_year(5:6,21); inc_year(4,21)],...%inc reduction in 2030
    [margins(1:2,[10:12,8],1); margins(5:6,[10:12,8],1); margins(4,[10:12,8],1)]]');%,...%prev in 2018, 2020, 2025
    %100*[R0_summary(1,1);R0_summary(1,2);0;0;R0_summary(1,7)]]'); 

paper_LB = round([total_treat_summary([1,2,5,6,3],:,2),...
    total_treat_2030_summary([1,2,5,6,3],4,2),...
    [treat_projected(3); sum(treat_projected(1:2));treat_scaleup(1,1);treat_scaleup(1,1);sum(treat_projected(1:2))],...
    [margins(1:2,1,2); margins(5:6,1,2); margins(3,1,2)],...
    [margins(1:2,2,2); margins(5:6,2,2); margins(3,2,2)],...
    [margins(1:2,[4,5,6,8,9],2); margins(5:6,[4,5,6,8,9],2); margins(3,[4,5,6,8,9],2)],...
    100-100*[inc_year_std(1:2,10,2); inc_year_std(5:6,10,2); inc_year_std(3,10,2)],...
    100-100*[inc_year_std(1:2,11,2); inc_year_std(5:6,11,2); inc_year_std(3,11,2)],...
    100-100*[inc_year_std(1:2,16,2); inc_year_std(5:6,16,2); inc_year_std(3,16,2)],...
    100-100*[inc_year_std(1:2,21,2); inc_year_std(5:6,21,2); inc_year_std(3,21,2)],...
    [margins(1:2,[10:12,8],2); margins(5:6,[10:12,8],2); margins(3,[10:12,8],2)]]');%,...
    %100*[R0_summary(2,1);R0_summary(2,2);0;0;R0_summary(2,7)]]');
paper_UB = round([total_treat_summary([1,2,5,6,3],:,3),...
    total_treat_2030_summary([1,2,5,6,3],4,3),...
    [treat_projected(3); sum(treat_projected(1:2));treat_scaleup(1,1);treat_scaleup(1,1);sum(treat_projected(1:2))],...
    [margins(1:2,1,3); margins(5:6,1,3); margins(3,1,3)],...
    [margins(1:2,2,3); margins(5:6,2,3); margins(3,2,3)],...
    [margins(1:2,[4,5,6,8,9],3); margins(5:6,[4,5,6,8,9],3); margins(3,[4,5,6,8,9],3)],...
    100-100*[inc_year_std(1:2,10,1); inc_year_std(5:6,10,1); inc_year_std(3,10,1)],...
    100-100*[inc_year_std(1:2,11,1); inc_year_std(5:6,11,1); inc_year_std(3,11,1)],...
    100-100*[inc_year_std(1:2,16,1); inc_year_std(5:6,16,1); inc_year_std(3,16,1)],...
    100-100*[inc_year_std(1:2,21,1); inc_year_std(5:6,21,1); inc_year_std(3,21,1)],...
    [margins(1:2,[10:12,8],3); margins(5:6,[10:12,8],3); margins(3,[10:12,8],3)]]');%,...
    %100*[R0_summary(3,1);R0_summary(3,2);0;0;R0_summary(3,7)]]');
paper_text = zeros(length(paper(:,1)),1+length(paper(1,:)));
paper_text = num2cell(paper_text);
paper_text{1,1} = 'Tr PWID';
paper_text{2,1} = 'Treat former PWID';
paper_text{3,1} = 'Treat other';
paper_text{4,1} = 'Treat total';
paper_text{5,1} = 'Treat total 2030';
paper_text{6,1} = 'Treat scale-up';
paper_text{7,1} = 'Cost';
paper_text{8,1} = 'QALYs';
paper_text{9,1} = 'ICER';
paper_text{10,1} = 'Incidence reduction';
paper_text{11,1} = 'Mortality reduction';
paper_text{12,1} = '2030 prevalence PWID';
paper_text{13,1} = 'Cumulative deaths';
paper_text{14,1} = 'Incidence reduction 2018';
paper_text{15,1} = 'Incidence reduction 2020';
paper_text{16,1} = 'Incidence reduction 2025';
paper_text{17,1} = 'Incidence reduction 2030';
paper_text{18,1} = 'PWID prevalence 2018';
paper_text{19,1} = 'PWID prevalence 2020';
paper_text{20,1} = 'PWID prevalence 2025';
paper_text{21,1} = 'PWID prevalence 2030';
paper_text{22,1} = 'R0';
for i = 1:length(paper(:,1))
    for j = 1:length(paper(1,:))
        paper_text{i,j+1} = [num2str(paper(i,j)),' (',num2str(paper_LB(i,j)),', ',num2str(paper_UB(i,j)),')'];
    end
end

