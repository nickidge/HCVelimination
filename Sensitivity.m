clear
global mu mu_PWID mu_former exit_IDU r_relapse delta alpha p_complete omega infect target_late prev0 age_cohort P...
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2...
    Q_sens q_svr q_treat sens c_daa discount total_PWID PWID0...
    c1_acutePWID c1_chronicPWID c2_PWID c3_PWID c4_PWID c5_PWID c6_PWID ...
    c1_acuteformerPWID c1_chronicformerPWID c2_formerPWID c3_formerPWID c4_formerPWID c5_formerPWID c6_formerPWID ...
    c1_acute c1_chronic c2 c3 c4 c5 c6 imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9...
    cascade0 cascade0_PWID disease0 cases0 cascade_target cascade_target_PWID scenario cascade_scale_time care age_mix...
    re_calibrate infected0 treat_no_old alpha_old Tin Run alpha_DAA harm_reduction treat n ...
    infect_base c1_chronicPWID_base c2_PWID_base c3_PWID_base c4_PWID_base c1_chronicformerPWID_base c2_formerPWID_base c3_formerPWID_base c4_formerPWID_base c1_chronic_base c2_base c3_base c4_base imp1_base imp2_base imp3_base imp4_base imp5_base imp6_base imp7_base imp8_base imp9_base ...
    target_inc target_death r_S4death_base start_year treat_projected

sens=1; %Number of runs in sensitivity analysis, sens=1 turns off feature
summary=zeros(6,12,sens);
R0_scens = zeros(sens,6);
Tin=66; % Years to run in model
Run=16; %Years to run the model
start_year = 20;
discount=0.03; % Discounting factor for costs and QALYs

mu=0.0083; %Average mortality rate for PWID
mu2024=-log(1-0.00044); %Mortality rates from Iceland life tables, assuming 66% male
mu2529=-log(1-0.00066);
mu3034=-log(1-0.00066);
mu3544=-log(1-0.00099);
mu4554=-log(1-0.00166);
mu5564=-log(1-0.00516);
mu6574=-log(1-0.01414);
mu7584=-log(1-0.04758);
mu85=-log(1-0.16045);
mu_former=[mu2024,mu2529,mu3034,mu3544,mu4554,mu5564,mu6574,mu7584,mu85];
mu2024=-log(1-0.0096); mu2529=-log(1-0.0096); mu3034=-log(1-0.0112); mu3545=-log(1-0.0018); %Mortality rates for PWID
mu_PWID=[mu2024,mu2529,mu3034,mu3544,mu4554,mu5564,mu6574,mu7584,mu85];

P=10; %Population to start the model
infected0=0.10; % initial proportion of PWID infected (in 1950)
PWID0=P*0.5; %Equilibrium proportion of PWID to former PWID
total_PWID = 400;

exit_IDU=1/14; %Cessation rate (1/length of injecting career)
r_relapse=-log(1-0.02); %Relapse to injecting rate

alpha_old=0.6; %SVR rate for interferon based therapies
alpha_DAA=0.90;% New treatment SVR rate, all genotypes
harm_reduction=0.0; %Harm reduction to reduce incidence that is also introduced

p_complete = 0.90; %Probability a PWID will complete treatment
omega=1/(18.36/52); % New treatment length, all genotypes
%omega=1/(24/52); % New treatment length, all genotypes
delta=0.26; %Proportion of infections that spontaneously clear

p_complete = 0.90; %Probability a PWID will complete treatment
omega=1/(18.36/52); % New treatment length, all genotypes
%omega=1/(24/52); % New treatment length, all genotypes
delta=0.26; %Proportion of infections that spontaneously clear

target_late=1; %Proportion of treatments allocated to late liver disease stage (set to -1 for proportionally)
treat_projected = [100,100,20]; %capped treatments; third entry is for old regimen
treat=treat_projected; 
age_cohort=23.5;
age_mix = 0.59; % proportion of injections occuring in same age bracket

imp1 = 0.001; % imported infections 1950-1975
imp2 = 0.001; % imported infections 1975-1980
imp3 = 0.001; % imported infections 1980-2085
imp4 = 5;
imp5 = 6;
imp6 = 15;
imp7 = 15;
imp8 = 7;
imp9 = 7;

target_inc=0.80; %Percentage reduction in incidence in WHO 2030 target
target_death=0.65; %Percentage reduction in deaths in WHO 2030 target

c_daa=[0,0]; %Cost of treatment is free (fixed total cost)

%% Calibration data
prev0 = xlsread('Template.xlsx','prev_PWID');
cascade0 = xlsread('Template.xlsx','cascade_all');
cascade0_PWID = xlsread('Template.xlsx','cascade_PWID');
disease0 = xlsread('Template.xlsx','disease');
cases0 = xlsread('Template.xlsx','cases'); % Antibody + diagnosed cases
cases0 = cases0(1:15,:);
HCC0 = xlsread('Template.xlsx','HCC'); % Antibody + diagnosed cases
diagnoses0 = xlsread('Template.xlsx','diagnoses'); % Antibody + diagnosed cases

%% Transition rates
r_AF0=52/12; %12 weeks
r_F0F1=-log(1-0.106);
r_F0F1_PWID=-log(1-0.116);
r_F1F2=-log(1-0.074);
r_F1F2_PWID=-log(1-0.085);
r_F2F3=-log(1-0.106);
r_F2F3_PWID=-log(1-0.085);
r_F3F4=-log(1-0.105);
r_F3F4_PWID=-log(1-0.13);
r_F4DC=-log(1-0.037);
r_DCHCC=-log(1-0.068);
r_F4HCC=-log(1-0.01);
r_HCCLT=-log(1-0.01);
r_DCLT=-log(1-0.0033);
r_DCdeath=-log(1-0.138);
r_HCCdeath=-log(1-0.605);
r_LTdeath1=-log(1-0.169);
r_LTdeath2=-log(1-0.034);
r_S4death=1/(-1/log(1-(0.138+0.605)/2)-1/log(1-0.01)); r_S4death_base = r_S4death; % 1% chance of HCC or DC following SVR at this stage
r_LT1LT2=1;

% Cascade of care parameters
c1_acutePWID = 1/20; % HCV antibody+ (acute phase)
c1_chronicPWID = 0.3; % HCV antibody+ (chronic phase)
c2_PWID = 0.4; % HCV RNA+
c3_PWID = 0.7; % Genotyped
c4_PWID = 0.9; % Liver disease tested
c5_PWID = 0.1; %Started treatment
c6_PWID = 1; % Not used yet
c1_acuteformerPWID = 0.1627; % HCV antibody+ (acute phase)
c1_chronicformerPWID = 0.65;
c2_formerPWID = 0.9;
c3_formerPWID = 0.73;
c4_formerPWID = 0.2;
c5_formerPWID = 0.1; %Started treatment
c6_formerPWID = 0.02; % Not used yet
c1_acute = 0.2099;
c1_chronic = 0.65;
c2 = 0.9;
c3 = 0.73;
c4 = 0.2;
c5 = 0.1; %Started treatment
c6 = 0.0; % Not used yet

cascade_scale_time = 0; % no scaling-up of cascade rates to start with
care = [0.3, 0.3]; % Proportion of specialist liver assessment at start and end of period
infect = 0.25; % initialize
re_calibrate=1;

%%
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point';
[summary_point]=sensitivity_func(filename)

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\annualtest_noRNA';
[summary_annualtest_noRNA]=sensitivity_func(filename)

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\2yeartest_noRNA';
[summary_2yeartest_noRNA]=sensitivity_func(filename)

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\6monthtest_noRNA';
[summary_6monthtest_noRNA]=sensitivity_func(filename)

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=1;
prev0 = [2015, 0.3; 2013, 0.3; 2010, 0.8; 2005, 0.3];
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_40prev';
[summary_40prev]=sensitivity_func(filename)
prev0 = xlsread('Template.xlsx','prev_PWID');

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=1;
prev0 = [2015, 0.5; 2013, 0.5; 2010, 0.5; 2005, 0.5];
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_60prev';
[summary_60prev]=sensitivity_func(filename)
prev0 = xlsread('Template.xlsx','prev_PWID')
re_calibrate = 0;

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=1;
total_PWID=200;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_60kPWID';
[summary_60kPWID]=sensitivity_func(filename)
total_PWID=400;

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=1;
total_PWID=1000;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_100kPWID';
[summary_100kPWID]=sensitivity_func(filename)
total_PWID=400;

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=1;
exit_IDU = 2 * exit_IDU;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_halfcareer';
[summary_halfcareer]=sensitivity_func(filename)
exit_IDU = 0.5 * exit_IDU;
re_calibrate=0;

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
%p_complete = 0.99; % Adherence among PWID
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\p_complete99';
[summary_p_comp99]=sensitivity_func(filename)

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\doubletime';
[summary_doubletime]=sensitivity_func(filename)

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
treat_projected = 2 * treat_projected;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_doubletreatnos';
[summary_doubletreatnos]=sensitivity_func(filename)
treat_projected = 0.5 * treat_projected;

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
alpha_DAA=0.95;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_99SVR';
[summary_99SVR]=sensitivity_func(filename)
alpha_DAA=0.90;

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
harm_reduction=0.0;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_0hr';
[summary_0hr]=sensitivity_func(filename)

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
harm_reduction=0.1;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_10hr';
[summary_10hr]=sensitivity_func(filename)

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
harm_reduction=0.2;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_20hr';
[summary_20hr]=sensitivity_func(filename)

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
harm_reduction=0.3;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_30hr';
[summary_30hr]=sensitivity_func(filename)

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=0;
harm_reduction=0.4;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_40hr';
[summary_40hr]=sensitivity_func(filename)

load('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
re_calibrate=1;
filename='C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_90diag';
[summary_90diag]=sensitivity_func(filename)
prev0 = xlsread('Template.xlsx','prev_PWID');
cascade0 = xlsread('Template.xlsx','cascade_all');
cascade0_PWID = xlsread('Template.xlsx','cascade_PWID');
disease0 = xlsread('Template.xlsx','disease');
cases0 = xlsread('Template.xlsx','cases'); % Antibody + diagnosed cases
cases0 = cases0(1:15,:);
HCC0 = xlsread('Template.xlsx','HCC'); % Antibody + diagnosed cases
diagnoses0 = xlsread('Template.xlsx','diagnoses'); % Antibody + diagnosed cases
