clear all
global filename scenario sens target_late Tin Run start_year num_pops num_cascade num_age num_intervention num_engagement num_region dt ... %meta parameters
    mu mu_PWID mu_former exit_IDU r_relapse delta alpha alpha_old alpha_DAA  omega age_cohort P y0 t0...
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2... % disease progression parameters
    Q_sens q_svr q_svr_PWID q_treat  c_daa discount treat  ... % QALY and cost parameters
    p_complete total_PWID PWID0 age_mix prev0 r_inc_up followup...%PWID parameters
    prev0_MSM ... %MSM parameters
    incarceration_rate prison_length ... %Prison parameters
    cascade0 cascade0_PWID cascade0_MSM disease0 disease0_HIV cases0 ost0 nsp0 HCC0 diagnoses0 infect_factor infect_base progression_base imported... % data and calibration
    ost_duration  nsp_duration treat_projected target_inc target_death cascade_scale_time care RNAtesting... %intervention
    infect progression imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9 ost_enrollment nsp_enrollment % calibtation parameters

user=extractBetween(pwd,"Users\","\");
drive=extractBefore(pwd,":");
filename=strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Testing\foo");

loaddata
dt = 1/4; % six-monthly time steps for burn-in / calibration perio 1950-2015
sens=1; %Number of runs in sensitivity analysis, sens=1 turns off feature
summary=zeros(6,12,sens);
followup = 1;


%% Sensitivity of parameters
s=1; % no sensitivity just yet

alpha=alpha_old; % calibrate in pre-DAA era
target_late=1; % target treatments to people with late-liver disease

%moved these up because need to use these for calibration if
%calibration is running.
infect_base=infect;
progression_base = progression;


[output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases,output_ost,output_nsp, output_diagnoses] = ...
    calibrate_optim_par(5, 4);
save('C:\Users\Nick\Desktop\Matlab Sims\Testing\prison_cal')
%load('calibration_data')
%filename is stored in calibration_data so have to add here so that we
%can have multiple users using these files
filename=strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Testing\foo");

infect_base=infect;
progression_base = progression;


%Run model in prior to 2016
[TT1,y1]=DE_track_age(Tin,y0,t0,treat);
y1_end=reshape(y1(end,:,:,:,:,:,:,:), num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region,33);
death2015 = (sum(sum(sum(sum(sum(sum(y1(end,:,:,:,:,:,:,22))))))) - sum(sum(sum(sum(sum(sum(y1(find(TT1>=65,1),:,:,:,:,:,:,22)))))))) / ...
    (TT1(end)-TT1(find(TT1>=65,1)));
inc2015 = sum(sum(sum(sum(sum(sum(sum(y1(find(TT1>=65,1):end-1,:,:,:,:,:,:,27))))))));%/(TT1(find(TT1>=66,1)-1) - TT1(find(TT1>=65,1)));
y1_end(:,:,:,:,:,:,21:(27+6))=0;
y1_end(:,6,:,:,:,:,1:20) = y1_end(:,6,:,:,:,:,1:20) + y1_end(:,8,:,:,:,:,1:20) + y1_end(:,10,:,:,:,:,1:20);
y1_end(:,8,:,:,:,:,1:20) = 0; y1_end(:,10,:,:,:,:,1:20) = 0; % moving failed to be re-eligible with DAAs


Prison_pops = reshape(sum(sum(sum(sum(sum(sum(y1(:,:,:,:,:,:,:,1:20),2),3),4),5),6),8),length(y1),num_region);
Prison_prev = reshape(sum(sum(sum(sum(sum(sum(y1(:,:,:,:,:,:,:,6:20),2),3),4),5),6),8),length(y1),num_region)./...
    reshape(sum(sum(sum(sum(sum(sum(y1(:,:,:,:,:,:,:,1:20),2),3),4),5),6),8),length(y1),num_region);

Prison_pops(:,2:end)
