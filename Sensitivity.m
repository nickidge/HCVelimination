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
directory=strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\");
calibration_file = strcat(directory,"calibration_draftv1.mat");


loaddata
dt = 1/4; % six-monthly time steps for burn-in / calibration perio 1950-2015
sens=1; %Number of runs in sensitivity analysis, sens=1 turns off feature
followup = 1;



%% Point
load(calibration_file)
filename=strcat(directory,'point');
[summary_point]=sens_func(filename)

%% double career
loaddata;
dt = 1/4; sens=1; summary=zeros(6,12,sens); followup = 1; alpha=alpha_old; target_late=1; infect_base=infect; progression_base = progression;
exit_IDU = 1/12;
[output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases,output_ost,output_nsp, output_diagnoses] = ...
    calibrate_optim_par(200, 30);
filename=strcat(directory,'HC');
save(strcat(directory,'HC_calibration'));
[summary_HC]=sens_func(filename)
exit_IDU = 1/6;

%% double prev
loaddata;
dt = 1/4; sens=1; summary=zeros(6,12,sens); followup = 1; alpha=alpha_old; target_late=1; infect_base=infect; progression_base = progression;
prev0(:,2) = 2 * prev0(:,2);
[output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases,output_ost,output_nsp, output_diagnoses] = ...
    calibrate_optim_par(200, 30);
filename=strcat(directory,'double_prev');
save(strcat(directory,'double_prev_calibration'));
[summary_double_prev]=sens_func(filename)

%% double PLHCV
loaddata;
dt = 1/4; sens=1; summary=zeros(6,12,sens); followup = 1; alpha=alpha_old; target_late=1; infect_base=infect; progression_base = progression;
cases0(:,2) = 1.2 * cases0(:,2);
[output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases,output_ost,output_nsp, output_diagnoses] = ...
    calibrate_optim_par(200, 30);
filename=strcat(directory,'double_PLHCV');
save(strcat(directory,'double_PLHCV_calibration'));
[summary_double_PLHCV]=sens_func(filename)

%% 70,000 PWID
loaddata;
dt = 1/4; sens=1; summary=zeros(6,12,sens); followup = 1; alpha=alpha_old; target_late=1; infect_base=infect; progression_base = progression;
total_PWID = 70,000;
[output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases,output_ost,output_nsp, output_diagnoses] = ...
    calibrate_optim_par(200, 30);
filename=strcat(directory,'70kPWID');
save(strcat(directory,'70kPWID_calibration'));
[summary_70kPWID]=sens_func(filename)

%% 30,000 PWID
loaddata;
dt = 1/4; sens=1; summary=zeros(6,12,sens); followup = 1; alpha=alpha_old; target_late=1; infect_base=infect; progression_base = progression;
total_PWID = 30,000;
[output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases,output_ost,output_nsp, output_diagnoses] = ...
    calibrate_optim_par(200, 30);
filename=strcat(directory,'30kPWID');
save(strcat(directory,'30kPWID_calibration'));
[summary_30kPWID]=sens_func(filename)