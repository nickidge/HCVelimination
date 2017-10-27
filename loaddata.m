function loaddata

global mu mu_PWID mu_former exit_IDU r_relapse delta alpha alpha_old alpha_DAA p_complete omega infect target_late dem prev0 prev0_MSM age_cohort P Tin Run y0_init y0 t0...
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2...
    Q_sens q_svr q_svr_PWID q_treat sens c_daa discount total_PWID PWID0 treat age_mix start_year...
    imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9 imported r_inc_up...
    cascade0 cascade0_PWID cascade0_MSM disease0 disease0_HIV cases0 ost0 nsp0 HCC0 diagnoses0 cascade_target cascade_target_PWID scenario cascade_scale_time care filename ...
    infect_base progression_base treat_projected ...
    num_pops num_cascade num_age num_intervention num_engagement num_region infect_factor progression...
    ost_enrollment ost_duration nsp_enrollment nsp_duration target_inc target_death 


%% Epi parameters
[epi0,epi0_label,epi0_raw] = xlsread('Template.xlsx','other_epi');
for i = 1:length(epi0_raw(2:end,1))
    str = epi0_raw(2+i-1,1);
    p = strcmpi(str,epi0_raw(:,1));
    epi0_raw(p,2);
    evalc(sprintf('%s=%s',cell2mat(str),num2str(cell2mat(epi0_raw(p,2)))));
end
alpha = alpha_old;
mu_former=[mu2024,mu2529,mu3034,mu3544,mu4554,mu5564,mu6574,mu7584,mu85];
mu_PWID=[mu2024_PWID,mu2529_PWID,mu3034_PWID,mu3544_PWID,mu4554,mu5564,mu6574,mu7584,mu85];
care = [care__a,care__b];
target_late=1; %Proportion of treatments allocated to late liver disease stage (set to -1 for proportionally)
treat_projected = [24000,16000,2000]; %capped treatments; third entry is for old regimen
treat=treat_projected;

dem = xlsread('Template.xlsx','demographics');

%% Model setup
[setup0,setup0_label,setup0_raw] = xlsread('Template.xlsx','structure');
for i = 1:length(setup0_raw(2:end,1))
    str = setup0_raw(2+i-1,1);
    p = strcmpi(str,setup0_raw(:,1));
    evalc(sprintf('%s=%s',cell2mat(str),num2str(cell2mat(setup0_raw(p,2)))));
end

infect_factor = zeros(num_pops,num_intervention, num_region);
%population group (current PWID, former PWID, other, MSM, MSM and current PWID, MSM and former PWID);
%intervention coverage (none, NSP, OST, NSP+OST);
infect_factor(:,:,1) = [[1,0,0];... 
    rel_incidence_NSP*[1,0,0];... %NSP
    rel_incidence_OST*[1,0,0];...%OST
    rel_incidence_NSPOST*[1,0,0]]'; %NSP + OST
%infect_factor(:,:,2) = infect_factor(:,:,1); %Assumnes prison the same but just no NSP
%infect_factor=1;  

%% Calibration data
prev0 = xlsread('Template.xlsx','prev_PWID');
%prev0_MSM = xlsread('Template.xlsx','prev_MSM');
cascade0 = xlsread('Template.xlsx','cascade_all');
cascade0_PWID = xlsread('Template.xlsx','cascade_PWID');
%cascade0_MSM = xlsread('Template.xlsx','cascade_MSM');
disease0 = xlsread('Template.xlsx','disease');
%disease0_HIV = xlsread('Template.xlsx','disease_HIV'); % Rate ratio for disease progression in HIV infected participants
cases0 = xlsread('Template.xlsx','cases'); % Antibody + diagnosed cases
[~,~,intervention_PWID0] = xlsread('Template.xlsx','Interventions_PWID'); % OST coverage
ost0 = [cell2mat(intervention_PWID0(2:end,1)), cell2mat(intervention_PWID0(2:end,strcmpi('OST',intervention_PWID0(1,:))))];
ost0 = ost0(~isnan(ost0(:,2)),:);
nsp0 = [cell2mat(intervention_PWID0(2:end,1)), cell2mat(intervention_PWID0(2:end,strcmpi('NSP',intervention_PWID0(1,:))))];
nsp0 = nsp0(~isnan(nsp0(:,2)),:);

HCC0 = xlsread('Template.xlsx','HCC'); % Antibody + diagnosed cases
diagnoses0 = xlsread('Template.xlsx','diagnoses'); % Antibody + diagnosed cases

%% Transition rates
[rates0,rates0_label,rates0_raw] = xlsread('Template.xlsx','transition_rates');

for i = 1:length(rates0_raw(2:end,1))
    str = rates0_raw(2+i-1,1);
    p = strcmpi(str,rates0_raw(:,1));
    evalc(sprintf('%s=%s',cell2mat(str),num2str(cell2mat(rates0_raw(p,2)))));
end


%% Calibration parameters initial guesses
[cal0,cal0_label,cal0_raw] = xlsread('Template.xlsx','calibration_guess');

for i = 1:length(cal0_raw(2:end,1))
    str = cal0_raw(2+i-1,1);
    p = strcmpi(str,cal0_raw(:,1));
    evalc(sprintf('%s=%s',cell2mat(str),num2str(cell2mat(cal0_raw(p,2)))));
end

% Cascade of care parameters
progression = zeros(num_pops,num_cascade,num_engagement,num_region);
progression(1,:,2,1) = [progression_PWID01,progression_PWID12,progression_PWID23,progression_PWID34,progression_PWID45,progression_PWID56,0,0,0,0]; %PWID
progression(2,:,2,1) = [progression_former01,progression_former12,progression_former23,progression_former34,progression_former45,progression_former56,0,0,0,0]; %Former PWID
progression(3,:,2,1) = [progression_other01,progression_other12,progression_other23,progression_other34,progression_other45,progression_other56,0,0,0,0]; %Other



%% Model setup
%Compartment dimensions for X(i,j,k,l,m):
%i, population group (current PWID, former PWID, other, MSM, MSM and current PWID, MSM and former PWID);
%j, cascade stage (undiagnosed, Ab diagnosed, RNA diagnosed, genotyped, liver assessed, on treatment, failed treatment, on re-treatment, failed re-treatment);
%k, age category
%l, intervention coverage status (none, NSP, OST, NSP+OST);
%m, care engagement (none, community-based, specialist)
%n, number of regions (community, prison)

S=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
S1=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
S2=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
S3=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
S4=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
A=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
T=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
T1=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
T2=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
T3=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
T4=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
F0=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
F1=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
F2=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
F3=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
F4=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
DC=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
HCC=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
LT=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
LT2=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
F4_transfer=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
Ldeath=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
T_total=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
HCC_transfer=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
T_F4on_total=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
Liver_transplants=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
Inc=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
Cas1=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
Cas2=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
Cas3=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
Cas4=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
Cas5=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);
Cas6=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region);

S(1,1,1,1,2,1)=(1-infected0)*PWID0;
S(2,1,1,1,2,1)=P-PWID0;
F0(1,1,1,1,2,1)=infected0*PWID0;


y0=reshape(cat(7,S,S1,S2,S3,S4,A,T,T1,T2,T3,T4,F0,F1,F2,F3,F4,DC,HCC,LT,LT2,F4_transfer,Ldeath,T_total,HCC_transfer,T_F4on_total,Liver_transplants,Inc,Cas1,Cas2,Cas3,Cas4,Cas5,Cas6),...
    num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region,33);
t0=0;

y0_init = y0;
end
