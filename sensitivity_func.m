
function[summary]=sensitivity_func(filename)

global mu mu_PWID mu_former exit_IDU r_relapse delta alpha p_complete omega infect target_late prev0 age_cohort P...
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2...
    Q_sens q_svr q_treat sens c_daa discount total_PWID PWID0...
    c1_acutePWID c1_chronicPWID c2_PWID c3_PWID c4_PWID c5_PWID c6_PWID ...
    c1_acuteformerPWID c1_chronicformerPWID c2_formerPWID c3_formerPWID c4_formerPWID c5_formerPWID c6_formerPWID ...
    c1_acute c1_chronic c2 c3 c4 c5 c6 imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9...
    cascade0 cascade0_PWID disease0 cases0 cascade_target cascade_target_PWID scenario cascade_scale_time care ...
    re_calibrate infected0 treat_no_old alpha_old Tin Run alpha_DAA harm_reduction treat n ...
    infect_base c1_chronicPWID_base c2_PWID_base c3_PWID_base c4_PWID_base c1_chronicformerPWID_base c2_formerPWID_base c3_formerPWID_base c4_formerPWID_base c1_chronic_base c2_base c3_base c4_base imp1_base imp2_base imp3_base imp4_base imp5_base imp6_base imp7_base imp8_base imp9_base ...
    target_inc target_death r_S4death_base treat_projected

summary=zeros(6,12,sens);

if strcmp(filename,'C:\Users\Nick\Desktop\Matlab Sims\Iceland\point_90diag')==1
    %% Calibration data
    prev0 = xlsread('Template_90diag.xlsx','prev_PWID');
    cascade0 = xlsread('Template_90diag.xlsx','cascade_all');
    cascade0_PWID = xlsread('Template_90diag.xlsx','cascade_PWID');
    disease0 = xlsread('Template_90diag.xlsx','disease');
    cases0 = xlsread('Template_90diag.xlsx','cases'); % Antibody + diagnosed cases
    cases0 = cases0(1:15,:);
    HCC0 = xlsread('Template_90diag.xlsx','HCC'); % Antibody + diagnosed cases
    diagnoses0 = xlsread('Template_90diag.xlsx','diagnoses'); % Antibody + diagnosed cases
end

%% Sensitivity of parameters
for s=1:sens
    s
    scenario = 'empty';
    cascade_scale_time = 0; % no scaling-up of cascade rates to start with
    care = [0.3, 0.3]; % Proportion of specialist liver assessment at start and end of period
    Q_sens=[0,0,0,0,0,0,0,0];
    q_svr=[0,0,0];
    q_treat=[0,0,0];
    if sens>1
        alpha=random('Uniform',0.9,0.9);
        p_complete=random('Uniform',0.85,0.95);
        
        r_AF0=52/random(truncate(makedist('Normal',12,2),1,26));
        r_F0F1=-log(1-random(truncate(makedist('Normal',0.106,.028),0.094,0.205)));
        r_F0F1_PWID=-log(1-random(truncate(makedist('Normal',0.116,.042),0.059,0.228)));
        r_F1F2=-log(1-random(truncate(makedist('Normal',0.074,.028),0.064,0.175)));
        r_F1F2_PWID=-log(1-random(truncate(makedist('Normal',0.085,.011),0.065,0.110)));
        r_F2F3=-log(1-random(truncate(makedist('Normal',0.106,.033),0.092,0.225)));
        r_F2F3_PWID=-log(1-random(truncate(makedist('Normal',0.085,.025),0.049,0.147)));
        r_F3F4=-log(1-random(truncate(makedist('Normal',0.105,.024),0.092,0.187)));
        r_F3F4_PWID=-log(1-random(truncate(makedist('Normal',0.13,.067),0.053,0.319)));
        r_F4DC=-log(1-random(truncate(makedist('Normal',0.037,.016),0.030,0.092)));
        r_DCHCC=-log(1-random(truncate(makedist('Normal',0.068,.015),0.041,0.099)));
        r_F4HCC=-log(1-random(truncate(makedist('Normal',0.01,.007),0.009,0.038)));
        r_HCCLT=-log(1-random(truncate(makedist('Normal',0.1,.033),0.050,0.180)));
        r_DCLT=-log(1-random(truncate(makedist('Normal',0.033,.008),0.017,0.049)));
        r_DCdeath=-log(1-random(truncate(makedist('Normal',0.138,.032),0.074,0.202)));
        r_HCCdeath=-log(1-random(truncate(makedist('Normal',0.605,.033),0.545,0.676)));
        r_LTdeath1=-log(1-random(truncate(makedist('Normal',0.169,.021),0.127,0.210)));
        r_LTdeath2=-log(1-random(truncate(makedist('Normal',0.034,.005),0.024,0.043)));
        r_S4death=-log(1-random(truncate(makedist('Normal',0.02,.001),0.01,0.03)));
        r_LT1LT2=1;
        q_S_PWID=random(truncate(makedist('Normal',0.930,.001),0.928,0.932)); %Have taken 0.1 off the S.D for q_A,F012,F3,F4,DC,HCC,LT1,LT2,svr(3),treat(1:3)
        q_S = 1;
        q_A=random(truncate(makedist('Normal',0.77,.012),0,q_S));
        q_F012=random(truncate(makedist('Normal',0.77,.012),0,q_S_PWID));
        q_F3=random(truncate(makedist('Normal',0.66,.015),0,q_F012));
        q_F4=random(truncate(makedist('Normal',0.55,.024),0,q_F3));
        q_DC=random(truncate(makedist('Normal',0.45,.014),0,q_F4));
        q_HCC=random(truncate(makedist('Normal',0.45,.014),0,q_F4));
        q_LT1=random(truncate(makedist('Normal',0.45,.014),q_HCC,1));
        q_LT2=random(truncate(makedist('Normal',0.67,.014),q_LT1,1));
        q_svr_PWID=[random(truncate(makedist('Normal',max(0.930,q_F012),.001),q_F012,q_S_PWID)),random(truncate(makedist('Normal',max(0.930,q_F012),.001),q_F012,q_S_PWID)),random(truncate(makedist('Normal',0.66,.05),q_DC,q_S_PWID))]; %After successful mild, standard and late treatment
        q_svr=[1,1,random(truncate(makedist('Normal',0.66,.05),q_DC,q_S))]; %After successful mild, standard and late treatment
        q_treat=[random(truncate(makedist('Normal',0.77,.012),q_F012,q_S)),random(truncate(makedist('Normal',0.66,.025),q_F3,q_S)),random(truncate(makedist('Normal',0.55,.024),q_DC,q_F3))]; %During mild, standard and late treatment
        delta=random('Uniform',0.22,0.29);
        Q_sens=[q_S_PWID,q_S,q_A,q_F012,q_F3,q_F4,q_DC,q_HCC,q_LT1,q_LT2];
    end
    %% Model setup
    %At equilibrium number of former PWID = mu2/(mu+r_relapse) * Current
    %Compartment dimensions represent X(i,j,k): i, injecting status (current, former, never); j, cascade stage; k, age category
    
    S(1,:,1)=[(1-infected0)*PWID0,zeros(1,9)];
    S(2,:,1)=[P-PWID0,zeros(1,9)];
    S(1:2,:,2:9)=zeros(2,10,8);
    S(3,:,:)=zeros(10,9);
    S1=zeros(3,10,9);
    S2=zeros(3,10,9);
    S3=zeros(3,10,9);
    S4=zeros(3,10,9);
    A=zeros(3,10,9);
    T=zeros(3,10,9);
    T1=zeros(3,10,9);
    T2=zeros(3,10,9);
    T3=zeros(3,10,9);
    T4=zeros(3,10,9);
    F0(1,:,1)=[infected0*PWID0,zeros(1,9)];
    F0(1,:,2:9)=zeros(10,8);
    F0(2:3,:,:)=zeros(2,10,9);
    F1=zeros(3,10,9);
    F2=zeros(3,10,9);
    F3=zeros(3,10,9);
    F4=zeros(3,10,9);
    DC=zeros(3,10,9);
    HCC=zeros(3,10,9);
    LT=zeros(3,10,9);
    LT2=zeros(3,10,9);
    F4_transfer=zeros(3,10,9);
    Ldeath=zeros(3,10,9);
    T_total=zeros(3,10,9);
    HCC_transfer=zeros(3,10,9);
    T_F4on_total=zeros(3,10,9);
    Liver_transplants=zeros(3,10,9);
    Inc=zeros(3,10,9);
    Cas1=zeros(3,10,9);
    Cas2=zeros(3,10,9);
    Cas3=zeros(3,10,9);
    Cas4=zeros(3,10,9);
    Cas5=zeros(3,10,9);
    Cas6=zeros(3,10,9);
    
    y0=cat(4,S,S1,S2,S3,S4,A,T,T1,T2,T3,T4,F0,F1,F2,F3,F4,DC,HCC,LT,LT2,F4_transfer,Ldeath,T_total,HCC_transfer,T_F4on_total,Liver_transplants,Inc,Cas1,Cas2,Cas3,Cas4,Cas5,Cas6);
    t0=0;
    %treat=[treat_no_old,0,0];
    
    alpha=alpha_old;
    target_late=1;
    %disease0 = 0;

    if re_calibrate==1
        [output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases] = calibrate_optim(40, 25);
        infect0=infect;
        c1_chronicPWID0 = c1_chronicPWID;
        c2_PWID0 = c2_PWID;
        c3_PWID0 = c3_PWID;
        c4_PWID0 = c4_PWID;
        c1_chronicformerPWID0 = c1_chronicformerPWID;
        c2_formerPWID0 = c2_formerPWID;
        c3_formerPWID0 = c3_formerPWID;
        c4_formerPWID0 = c4_formerPWID;
        c1_chronic0 = c1_chronic;
        c20 = c2;
        c30 = c3;
        c40 = c4;
        imp10 = imp1;
        imp20 = imp2;
        imp30 = imp3;
        imp40 = imp4;
        imp50 = imp5;
        imp60 = imp6;
        imp70 = imp7;
        imp80 = imp8;
        imp90 = imp9;
        r_S4death0 = r_S4death_base;
        if strcmp(filename,'C:\Users\Nick\Desktop\Matlab Sims\Iceland\point')==1
            
            infect_base=infect0;
            c1_chronicPWID_base = c1_chronicPWID0;
            c2_PWID_base = c2_PWID0;
            c3_PWID_base = c3_PWID0;
            c4_PWID_base = c4_PWID0;
            c1_chronicformerPWID_base = c1_chronicformerPWID0;
            c2_formerPWID_base = c2_formerPWID0;
            c3_formerPWID_base = c3_formerPWID0;
            c4_formerPWID_base = c4_formerPWID0;
            c1_chronic_base = c1_chronic0;
            c2_base = c20;
            c3_base = c30;
            c4_base = c40;
            imp1_base = imp10;
            imp2_base = imp20;
            imp3_base = imp30;
            imp4_base = imp40;
            imp5_base = imp50;
            imp6_base = imp60;
            imp7_base = imp70;
            imp8_base = imp80;
            imp9_base = imp90;
            %r_S4death_base = r_S4death0;
            save('C:\Users\Nick\Desktop\Matlab Sims\Iceland\calibrate');
        end
    else
        infect0=infect_base;
        c1_chronicPWID0 = c1_chronicPWID_base;
        c2_PWID0 = c2_PWID_base;
        c3_PWID0 = c3_PWID_base;
        c4_PWID0 = c4_PWID_base;
        c1_chronicformerPWID0 = c1_chronicformerPWID_base;
        c2_formerPWID0 = c2_formerPWID_base;
        c3_formerPWID0 = c3_formerPWID_base;
        c4_formerPWID0 = c4_formerPWID_base;
        c1_chronic0 = c1_chronic_base;
        c20 = c2_base;
        c30 = c3_base;
        c40 = c4_base;
        imp10 = imp1_base;
        imp20 = imp2_base;
        imp30 = imp3_base;
        imp40 = imp4_base;
        imp50 = imp5_base;
        imp60 = imp6_base;
        imp70 = imp7_base;
        imp80 = imp8_base;
        imp90 = imp9_base;
        r_S4death0 = r_S4death_base;
        infect=infect_base;
        c1_chronicPWID = c1_chronicPWID_base;
        c2_PWID = c2_PWID_base;
        c3_PWID = c3_PWID_base;
        c4_PWID = c4_PWID_base;
        c1_chronicformerPWID = c1_chronicformerPWID_base;
        c2_formerPWID = c2_formerPWID_base;
        c3_formerPWID = c3_formerPWID_base;
        c4_formerPWID = c4_formerPWID_base;
        c1_chronic = c1_chronic_base;
        c2 = c2_base;
        c3 = c3_base;
        c4 = c4_base;
        imp1 = imp1_base;
        imp2 = imp2_base;
        imp3 = imp3_base;
        imp4 = imp4_base;
        imp5 = imp5_base;
        imp6 = imp6_base;
        imp7 = imp7_base;
        imp8 = imp8_base;
        imp9 = imp9_base;
        r_S4death = r_S4death_base;
    end
    
    %Run model in prior to 2016
    [TT1,y1]=DE_track_age(Tin,y0,t0,treat);
    y1_end=reshape(y1(end,:,:,:,:),3,10,9,27+6);
    death2015 = (sum(sum(sum(y1(end,:,:,:,22)))) - sum(sum(sum(y1(find(TT1>=65,1),:,:,:,22))))) / ...
        (TT1(end)-TT1(find(TT1>=65,1)));
    inc2015 = sum(sum(sum(sum(y1(find(TT1>=65,1)+1:end,:,:,:,27)))))/(66 - TT1(find(TT1>=65,1)+1));
    y1_end(:,:,:,21:(27+6))=0;
    y1_end(:,6,:,1:20) = y1_end(:,6,:,1:20) + y1_end(:,8,:,1:20) + y1_end(:,10,:,1:20); y1_end(:,8,:,1:20) = 0; y1_end(:,10,:,1:20) = 0; % moving failed to be re-eligible with DAAs
    
    if strcmp(filename,'C:\Users\Nick\Desktop\Matlab Sims\Iceland\rapidRNA')==1
        c2_PWID = 100; c2_formerPWID = 100; c2 = 100;
    end
    if strcmp(filename,'C:\Users\Nick\Desktop\Matlab Sims\Iceland\r_S4death005')==1
        r_S4death = 0.005;
    end
    if strcmp(filename,'C:\Users\Nick\Desktop\Matlab Sims\Iceland\p_complete99')==1
        p_complete = 0.99;
    end
    
    %% Baseline: Current standard of care with no scaled up treatment
    scenario = 'base'; %Current level of community care
    c5_PWID = 1/(30/365); c5_formerPWID = 1/(30/365); c5 = 1/(30/365);
    c6_PWID = 1/(30/365); c6_formerPWID = 1/(30/365); c6 = 1/(30/365);
    treat=[20,0,20]; %capped treatments; third entry is for old regimen
    target_late=1; alpha = alpha_DAA;
    infect = infect * (1-harm_reduction);
    [TT2,y2]=DE_track_age(Run,y1_end,TT1,treat);
    
    ycomb=[y1(1:end-1,:,:,:,:);y2];
    %Make everyone susceptible for R0 calculation
    ycomb_S = zeros(length(TT2),3,10,9,20);
    ycomb_S(:,:,1,:,1) = ycomb(:,:,1,:,1) + sum(sum(ycomb(:,:,2:10,:,[6,7,12]),3),5);
    ycomb_S(:,:,1,:,2) = ycomb(:,:,1,:,2) + sum(sum(ycomb(:,:,2:10,:,[8,13]),3),5);
    ycomb_S(:,:,1,:,3) = ycomb(:,:,1,:,3) + sum(sum(ycomb(:,:,2:10,:,[9,14]),3),5);
    ycomb_S(:,:,1,:,4) = ycomb(:,:,1,:,4) + sum(sum(ycomb(:,:,2:10,:,[10,15]),3),5);
    ycomb_S(:,:,1,:,5) = ycomb(:,:,1,:,5) + sum(sum(ycomb(:,:,2:10,:,[11,16:20]),3),5);
    R0_scens(1) = R0(65,real(reshape(ycomb_S(find(TT2>65,1),:,:,:,1:20),3*10*9*20,1)))
    
    y2_noage=reshape(sum(y2(:,:,:,:,:),4),length(y2(:,1,1,1,1)),3,10,27+6); %Reshape to sum over age stratification
    ycomb_noage=reshape(sum(ycomb(:,:,:,:,:),4),length(ycomb(:,1,1,1,1)),3,10,27+6); %Reshape to sum over age stratification
    
    Tint=Tin+[0,3,14]; t_val=zeros(1,length(Tint));
    for l=1:length(Tint) %Find times corresponding to years since intervention
        t_val(l)=find(TT2>=(Tint(l)),1);
    end
    TT_=TT2(t_val(1):t_val(2))-Tin;
    [cost,qaly,life_exp,tr]=Costs_age(y2_noage(1:t_val(2)-t_val(1)+1,:,:,:),TT_);
    tr_(1) = sum(sum(sum(y2_noage(1:t_val(3)-t_val(1)+1,1,:,[23,25])))); tr_(2) = sum(sum(sum(y2_noage(1:t_val(3)-t_val(1)+1,1,:,[23,25])))); tr_(3) = sum(sum(sum(y2_noage(1:t_val(3)-t_val(1)+1,3,:,[23,25])))); tr_(4) = sum(sum(sum(sum(y2_noage(1:t_val(3)-t_val(1)+1,:,:,[23,25])))));
    summary(1,:,s)=[cost,qaly,tr(4),tr_(4),...
        sum(sum(sum(ycomb_noage(find(TT2>=80,1)+1:find(TT2>=81,1),:,:,27))))/(TT2(find(TT2>=81,1))-TT2(find(TT2>=80,1))),... %Incidence in 2030
        (sum(sum(ycomb_noage(find(TT2>=81,1),1:3,:,22)))-sum(sum(ycomb_noage(find(TT2>=80,1),1:3,:,22))))...
        /(TT2(find(TT2>=81,1))-TT2(find(TT2>=80,1))),... %Liver related deaths in 2030
        100*sum(sum(ycomb_noage(t_val(1),1,:,[6,12:20])))./sum(sum(ycomb_noage(t_val(1),1,:,1:20))),... %Prevalence among PWID in 2015
        100*sum(sum(ycomb_noage(t_val(3),1,:,[6,12:20])))./sum(sum(ycomb_noage(t_val(3),1,:,1:20))),... %Prevalence among PWID in 2030
        sum(sum(ycomb_noage(t_val(3),1:3,:,22))),... %Liver related deaths 2015 to 2030
        100*sum(sum(ycomb_noage(find(TT2>=69,1),1,:,[6,12:20])))./sum(sum(ycomb_noage(find(TT2>=69,1),1,:,1:20))),... %Prevalence among PWID in 2019
        100*sum(sum(ycomb_noage(find(TT2>=70,1),1,:,[6,12:20])))./sum(sum(ycomb_noage(find(TT2>=70,1),1,:,1:20))),... %Prevalence among PWID in 2020
        100*sum(sum(ycomb_noage(find(TT2>=75,1),1,:,[6,12:20])))./sum(sum(ycomb_noage(find(TT2>=75,1),1,:,1:20)))]; %Prevalence among PWID in 2025
       
    
    %%  Scenario 1: Scale up treatment as projected
    scenario = 'current'
    c5_PWID = 1/(30/365); c5_formerPWID = 1/(30/365); c5 = 1/(30/365);
    c6_PWID = 1/(30/365); c6_formerPWID = 1/(30/365); c6 = 1/(30/365);
    if strcmp(filename,'C:\Users\Nick\Desktop\Matlab Sims\Iceland\doubletime')==1
        c5_PWID = 1/(2*30/365); c5_formerPWID = 1/(2*30/365); c5 = 1/(2*30/365);
        c6_PWID = 1/(2*30/365); c6_formerPWID = 1/(2*30/365); c6 = 1/(2*30/365);
        c5_PWID0 = 1/(2*30/365); c5_formerPWID0 = 1/(2*30/365); c50 = 1/(2*30/365);
        c6_PWID0 = 1/(2*30/365); c6_formerPWID0 = 1/(2*30/365); c60 = 1/(2*30/365);
    end
    if strcmp(filename,'C:\Users\Nick\Desktop\Matlab Sims\Iceland\2yeartest_noRNA')==1
        temp_rate = (c1_chronicPWID_base - 0.9*0.5) / (0.1);
        c1_chronicPWID = 0.5;
        c1_chronicPWID0 = 0.5;
        c1_chronicformerPWID = 0.5;
        c1_chronicformerPWID0 = 0.5;
    end
    if strcmp(filename,'C:\Users\Nick\Desktop\Matlab Sims\Iceland\6monthtest_noRNA')==1
        temp_rate = (c1_chronicPWID_base - 0.9*0.5) / (0.1);
        c1_chronicPWID = 2;
        c1_chronicPWID0 = 2;
        c1_chronicformerPWID = 2;
        c1_chronicformerPWID0 = 2;
    end
    if strcmp(filename,'C:\Users\Nick\Desktop\Matlab Sims\Iceland\annualtest_noRNA')==1
        temp_rate = (c1_chronicPWID - 0.31*1) / (1-0.31);
        c1_chronicPWID = 1;
        c1_chronicPWID0 = 1;
        c1_chronicformerPWID = 1;
        c1_chronicformerPWID0 = 1;
    end
    target_late=-1; alpha = alpha_DAA; treat=treat_projected;
    [TT2_treat,y2_treat]=DE_track_age(Run,y1_end,TT1,treat);
    
    ycomb2=[y1(1:end-1,:,:,:,:);y2_treat];
    y2_treat_noage=reshape(sum(y2_treat(:,:,:,:,:),4),length(y2_treat(:,1,1,1,1)),3,10,27+6);
    ycomb2_noage=reshape(sum(ycomb2(:,:,:,:,:),4),length(ycomb2(:,1,1,1,1)),3,10,27+6);
    
    ycomb2_S = zeros(length(TT2_treat),3,10,9,20);
    ycomb2_S(:,:,1,:,1) = ycomb2(:,:,1,:,1) + sum(sum(ycomb2(:,:,2:10,:,[6,7,12]),3),5);
    ycomb2_S(:,:,1,:,2) = ycomb2(:,:,1,:,2) + sum(sum(ycomb2(:,:,2:10,:,[8,13]),3),5);
    ycomb2_S(:,:,1,:,3) = ycomb2(:,:,1,:,3) + sum(sum(ycomb2(:,:,2:10,:,[9,14]),3),5);
    ycomb2_S(:,:,1,:,4) = ycomb2(:,:,1,:,4) + sum(sum(ycomb2(:,:,2:10,:,[10,15]),3),5);
    ycomb2_S(:,:,1,:,5) = ycomb2(:,:,1,:,5) + sum(sum(ycomb2(:,:,2:10,:,[11,16:20]),3),5);
    R0_scens(2) = R0(68,real(reshape(ycomb2_S(find(TT2_treat>68,1),:,:,:,1:20),3*10*9*20,1)));
    
    Tint=Tin+[0,3,14]; t_val_treat=zeros(1,length(Tint));
    for l=1:length(Tint) %Find times corresponding to years since intervention
        t_val_treat(l)=find(TT2_treat>=(Tint(l)),1);
    end
    TT_2=TT2_treat(t_val_treat(1):t_val_treat(2))-Tin;
    
    
    [cost_treat,qaly_treat,life_exp_treat,tr2]=Costs_age(y2_treat_noage(1:t_val_treat(2)-t_val_treat(1)+1,:,:,:),TT_2);
    tr2_(1) = sum(sum(sum(y2_treat_noage(1:t_val_treat(3)-t_val_treat(1)+1,1,:,[23,25])))); tr2_(2) = sum(sum(sum(y2_treat_noage(1:t_val_treat(3)-t_val_treat(1)+1,1,:,[23,25])))); tr2_(3) = sum(sum(sum(y2_treat_noage(1:t_val_treat(3)-t_val_treat(1)+1,3,:,[23,25])))); tr2_(4) = sum(sum(sum(sum(y2_treat_noage(1:t_val_treat(3)-t_val_treat(1)+1,:,:,[23,25])))));
    summary(2,:,s)=[cost_treat,qaly_treat,tr2(4),tr2_(4),...
        sum(sum(sum(ycomb2_noage(find(TT2_treat>=80,1)+1:find(TT2_treat>=81,1),:,:,27))))/(TT2_treat(find(TT2_treat>=81,1))-TT2_treat(find(TT2_treat>=80,1))),... %Incidence in 2030
        (sum(sum(ycomb2_noage(find(TT2_treat>=81,1),1:3,:,22)))-sum(sum(ycomb2_noage(find(TT2_treat>=80,1),1:3,:,22))))...
        /(TT2_treat(find(TT2_treat>=81,1))-TT2_treat(find(TT2_treat>=80,1))),... %Liver related deaths in 2030
        100*sum(sum(ycomb2_noage(t_val_treat(1)-1,1,:,[6,12:20])))./sum(sum(ycomb2_noage(t_val_treat(1),1,:,1:20))),... %Prevalence among PWID in 2015
        100*sum(sum(ycomb2_noage(t_val_treat(3),1,:,[6,12:20])))./sum(sum(ycomb2_noage(t_val_treat(3),1,:,1:20))),...%Prevalence among PWID in 2030
        sum(sum(ycomb2_noage(t_val_treat(3),1:3,:,22))),... %Liver related deaths 2015 to 2030
        100*sum(sum(ycomb2_noage(find(TT2_treat>=69,1),1,:,[6,12:20])))./sum(sum(ycomb2_noage(find(TT2_treat>=69,1),1,:,1:20))),... %Prevalence among PWID in 2019
        100*sum(sum(ycomb2_noage(find(TT2_treat>=70,1),1,:,[6,12:20])))./sum(sum(ycomb2_noage(find(TT2_treat>=70,1),1,:,1:20))),... %Prevalence among PWID in 2020
        100*sum(sum(ycomb2_noage(find(TT2_treat>=75,1),1,:,[6,12:20])))./sum(sum(ycomb2_noage(find(TT2_treat>=75,1),1,:,1:20)))]; %Prevalence among PWID in 2025
    
    %%  Scenario 2: Calculate required treatments for incidence target    
    scenario = 'current'
    target_late=0; % Target PWID
    alpha = alpha_DAA;
    %apply_cascade_rates(scenario, 0);
    c1_chronicPWID = 50; c1_chronicformerPWID = 50; c1_chronic = 50;
    c2_PWID = 50; c2_formerPWID = 50; c2 = 50;
    c3_PWID = 50; c3_formerPWID = 50; c3 = 50;
    c4_PWID = 50; c4_formerPWID = 50; c4 = 50;
    c5_PWID = 50; c5_formerPWID = 50; c5 = 50;
    c6_PWID = 50; c6_formerPWID = 50; c6 = 50;
    %scenario = 'rapidRNA'; c2_PWID = 100; c2_formerPWID = 100; c2 = 100;
    %infect = (1-0.25) * infect/(1-harm_reduction); % Harm reduction at 20% scaleup instead
    
    [TT2_treat5,y2_treat5]=DE_track_age(Run,y1_end,TT1,treat);
    elim_years = [69,70,75,80];
    if strcmp(scenario,'current') == 1
        for j = 1:length(elim_years)
            treat_range = 5:5:565; % Estimate of testing required
            for i = 1:length(treat_range)
                i
                treat(1) = treat_range(i); treat(2) = 0;
                [TT2_treat5,y2_treat5]=DE_track_age(Run,y1_end,TT1,treat);
                ycomb5=[y1(1:end-1,:,:,:,:);y2_treat5];
                y2_treat5_noage=reshape(sum(y2_treat5(:,:,:,:,:),4),length(y2_treat5(:,1,1,1,1)),3,10,27+6);
                ycomb5_noage=reshape(sum(ycomb5(:,:,:,:,:),4),length(ycomb5(:,1,1,1,1)),3,10,27+6);
                incred(i,j) = (inc2015 - sum(sum(sum(ycomb5_noage(find(TT2_treat5>=elim_years(j),1)+1:find(TT2_treat5>=elim_years(j)+1,1),:,:,27))))/...
                    (TT2_treat5(find(TT2_treat5>=elim_years(j)+1,1))-TT2_treat5(find(TT2_treat5>=elim_years(j),1)))) / inc2015;
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
    ycomb5=[y1(1:end-1,:,:,:,:);y2_treat5];
    y2_treat5_noage=reshape(sum(y2_treat5(:,:,:,:,:),4),length(y2_treat5(:,1,1,1,1)),3,10,27+6);
    ycomb5_noage=reshape(sum(ycomb5(:,:,:,:,:),4),length(ycomb5(:,1,1,1,1)),3,10,27+6);
    
    Tint=Tin+[0,3,14]; t_val_treat5=zeros(1,length(Tint));
    for l=1:length(Tint) %Find times corresponding to years since intervention
        t_val_treat5(l)=find(TT2_treat5>=(Tint(l)),1);
    end
    TT_5=TT2_treat5(t_val_treat5(1):t_val_treat5(2))-Tin;
    
    
    [cost_treat5,qaly_treat5,life_exp_treat5,tr5]=Costs_age(y2_treat5_noage(1:t_val_treat5(2)-t_val_treat5(1)+1,:,:,:),TT_5);
    tr5_(1) = sum(sum(sum(y2_treat5_noage(1:t_val_treat5(3)-t_val_treat5(1)+1,1,:,[23,25])))); tr5_(2) = sum(sum(sum(y2_treat5_noage(1:t_val_treat5(3)-t_val_treat5(1)+1,1,:,[23,25])))); tr5_(3) = sum(sum(sum(y2_treat5_noage(1:t_val_treat5(3)-t_val_treat5(1)+1,3,:,[23,25])))); tr5_(4) = sum(sum(sum(sum(y2_treat5_noage(1:t_val_treat5(3)-t_val_treat5(1)+1,:,:,[23,25])))));
    summary(5,:,s)=[cost_treat5,qaly_treat5,tr5(4),tr5_(4),...
        sum(sum(sum(ycomb5_noage(find(TT2_treat5>=80,1)+1:find(TT2_treat5>=81,1),:,:,27))))/(TT2_treat5(find(TT2_treat5>=81,1))-TT2_treat5(find(TT2_treat5>=80,1))),... %Incidence in 2030
        (sum(sum(ycomb5_noage(find(TT2_treat5>=81,1),1:3,:,22)))-sum(sum(ycomb5_noage(find(TT2_treat5>=80,1),1:3,:,22))))...
        /(TT2_treat5(find(TT2_treat5>=81,1))-TT2_treat5(find(TT2_treat5>=80,1))),... %Liver related deaths in 2030
        100*sum(sum(ycomb5_noage(t_val_treat5(1)-1,1,:,[6,12:20])))./sum(sum(ycomb5_noage(t_val_treat5(1),1,:,1:20))),... %Prevalence among PWID in 2015
        100*sum(sum(ycomb5_noage(t_val_treat5(3),1,:,[6,12:20])))./sum(sum(ycomb5_noage(t_val_treat5(3),1,:,1:20))),...%Prevalence among PWID in 2030
        sum(sum(ycomb5_noage(t_val_treat5(3),1:3,:,22))),... %Liver related deaths 2015 to 2030
        100*sum(sum(ycomb5_noage(find(TT2_treat5>=69,1),1,:,[6,12:20])))./sum(sum(ycomb5_noage(find(TT2_treat5>=69,1),1,:,1:20))),... %Prevalence among PWID in 2019
        100*sum(sum(ycomb5_noage(find(TT2_treat5>=70,1),1,:,[6,12:20])))./sum(sum(ycomb5_noage(find(TT2_treat5>=70,1),1,:,1:20))),... %Prevalence among PWID in 2020
        100*sum(sum(ycomb5_noage(find(TT2_treat5>=75,1),1,:,[6,12:20])))./sum(sum(ycomb5_noage(find(TT2_treat5>=75,1),1,:,1:20)))]; %Prevalence among PWID in 2025

    %%  Scenario 3: Calculate required treatments for mortality target
    scenario = 'current'
    target_late=0; % Target PWID
    alpha = alpha_DAA;
    c1_chronicPWID = 50; c1_chronicformerPWID = 50; c1_chronic = 50;
    c2_PWID = 50; c2_formerPWID = 50; c2 = 50;
    c3_PWID = 50; c3_formerPWID = 50; c3 = 50;
    c4_PWID = 50; c4_formerPWID = 50; c4 = 50;
    c5_PWID = 50; c5_formerPWID = 50; c5 = 50;
    c6_PWID = 50; c6_formerPWID = 50; c6 = 50;
    
    treat(1) = treat_scaleup(2,s); treat(2) = 0;
    
    [TT2_treat6,y2_treat6]=DE_track_age(Run,y1_end,TT1,treat);
    ycomb6=[y1(1:end-1,:,:,:,:);y2_treat6];
    y2_treat6_noage=reshape(sum(y2_treat6(:,:,:,:,:),4),length(y2_treat6(:,1,1,1,1)),3,10,27+6);
    ycomb6_noage=reshape(sum(ycomb6(:,:,:,:,:),4),length(ycomb6(:,1,1,1,1)),3,10,27+6);
    
    Tint=Tin+[0,3,14]; t_val_treat6=zeros(1,length(Tint));
    for l=1:length(Tint) %Find times corresponding to years since intervention
        t_val_treat6(l)=find(TT2_treat6>=(Tint(l)),1);
    end
    TT_6=TT2_treat6(t_val_treat6(1):t_val_treat6(2))-Tin;
    
    
    [cost_treat6,qaly_treat6,life_exp_treat6,tr6]=Costs_age(y2_treat6_noage(1:t_val_treat6(2)-t_val_treat6(1)+1,:,:,:),TT_6);
    tr6_(1) = sum(sum(sum(y2_treat6_noage(1:t_val_treat6(3)-t_val_treat6(1)+1,1,:,[23,25])))); tr6_(2) = sum(sum(sum(y2_treat6_noage(1:t_val_treat6(3)-t_val_treat6(1)+1,1,:,[23,25])))); tr6_(3) = sum(sum(sum(y2_treat6_noage(1:t_val_treat6(3)-t_val_treat6(1)+1,3,:,[23,25])))); tr6_(4) = sum(sum(sum(sum(y2_treat6_noage(1:t_val_treat6(3)-t_val_treat6(1)+1,:,:,[23,25])))));
    summary(6,:,s)=[cost_treat6,qaly_treat6,tr6(4),tr6_(4),...
        sum(sum(sum(ycomb6_noage(find(TT2_treat6>=80,1)+1:find(TT2_treat6>=81,1),:,:,27))))/(TT2_treat6(find(TT2_treat6>=81,1))-TT2_treat6(find(TT2_treat6>=80,1))),... %Incidence in 2030
        (sum(sum(ycomb6_noage(find(TT2_treat6>=81,1),1:3,:,22)))-sum(sum(ycomb6_noage(find(TT2_treat6>=80,1),1:3,:,22))))...
        /(TT2_treat6(find(TT2_treat6>=81,1))-TT2_treat6(find(TT2_treat6>=80,1))),... %Liver related deaths in 2030
        100*sum(sum(ycomb6_noage(t_val_treat6(1)-1,1,:,[6,12:20])))./sum(sum(ycomb6_noage(t_val_treat6(1),1,:,1:20))),... %Prevalence among PWID in 2015
        100*sum(sum(ycomb6_noage(t_val_treat6(3),1,:,[6,12:20])))./sum(sum(ycomb6_noage(t_val_treat6(3),1,:,1:20))),...%Prevalence among PWID in 2030
        sum(sum(ycomb6_noage(t_val_treat6(3),1:3,:,22))),... %Liver related deaths 2015 to 2030
        100*sum(sum(ycomb6_noage(find(TT2_treat6>=69,1),1,:,[6,12:20])))./sum(sum(ycomb6_noage(find(TT2_treat6>=69,1),1,:,1:20))),... %Prevalence among PWID in 2019
        100*sum(sum(ycomb6_noage(find(TT2_treat6>=70,1),1,:,[6,12:20])))./sum(sum(ycomb6_noage(find(TT2_treat6>=70,1),1,:,1:20))),... %Prevalence among PWID in 2020
        100*sum(sum(ycomb6_noage(find(TT2_treat6>=75,1),1,:,[6,12:20])))./sum(sum(ycomb6_noage(find(TT2_treat6>=75,1),1,:,1:20)))]; %Prevalence among PWID in 2025
    
    %%  Scenario 4: Testing required to get R0 < 1
    scenario = 'current';
    c1_chronicPWID = c1_chronicPWID0; c1_chronicformerPWID = c1_chronicformerPWID0; c1_chronic = c1_chronic0;
    c2_PWID = c2_PWID0; c2_formerPWID = c2_formerPWID0; c2 = c20;
    c3_PWID = c3_PWID0; c3_formerPWID = c3_formerPWID0; c3 = c30;
    c4_PWID = c4_PWID0; c4_formerPWID = c4_formerPWID0; c4 = c40;
    c5_PWID = 1/(30/365); c5_formerPWID = 1/(30/365); c5 = 1/(30/365);
    c6_PWID = 1/(30/365); c6_formerPWID = 1/(30/365); c6 = 1/(30/365);
    if strcmp(filename,'C:\Users\Nick\Desktop\Matlab Sims\Iceland\doubletime')==1
        c5_PWID = 1/(2*30/365); c5_formerPWID = 1/(2*30/365); c5 = 1/(2*30/365);
        c6_PWID = 1/(2*30/365); c6_formerPWID = 1/(2*30/365); c6 = 1/(2*30/365);
    end
    target_late = -1;
    range_test = [-1, 1/2, 2, 1];
    infect_red = 0;
    treat = treat_projected;
    for i = 1:length(range_test)
        for j = 1:length(infect_red)
            if range_test(i) > 0
                temp_rate = (c1_chronicPWID_base - 0.9*0.5) / (0.1);
                c1_chronicPWID = 0.9*range_test(i) + 0.1*temp_rate;
                c1_chronicformerPWID = 0.9*range_test(i) + 0.1*temp_rate;
                %c1_chronicPWID = range_test(i);
                %c1_chronicformerPWID = range_test(i);
            end
            infect = infect_base * (1-infect_red(j));
            %if c2_PWID < 100 c2_PWID = range(i) * c2_PWID_base; end
            
            [TT2_treat3,y2_treat3]=DE_track_age(Run,y1_end,TT1,treat);
            ycomb3=[y1(1:end-1,:,:,:,:);y2_treat3];
            
            y2_treat3_noage=reshape(sum(y2_treat3(:,:,:,:,:),4),length(y2_treat3(:,1,1,1,1)),3,10,27+6);
            ycomb3_noage=reshape(sum(ycomb3(:,:,:,:,:),4),length(ycomb3(:,1,1,1,1)),3,10,27+6);
            Tint=Tin+[0,3,14]; t_val_treat3=zeros(1,length(Tint));
            for l=1:length(Tint) %Find times corresponding to years since intervention
                t_val_treat3(l)=find(TT2_treat3>=(Tint(l)),1);
            end
            TT_3=TT2_treat3(t_val_treat3(1):t_val_treat3(2))-Tin;
            
            ycomb3_S = zeros(length(TT2_treat3),3,10,9,20);
            ycomb3_S(:,:,1,:,1) = ycomb3(:,:,1,:,1) + sum(sum(ycomb3(:,:,2:10,:,[6,7,12]),3),5);
            ycomb3_S(:,:,1,:,2) = ycomb3(:,:,1,:,2) + sum(sum(ycomb3(:,:,2:10,:,[8,13]),3),5);
            ycomb3_S(:,:,1,:,3) = ycomb3(:,:,1,:,3) + sum(sum(ycomb3(:,:,2:10,:,[9,14]),3),5);
            ycomb3_S(:,:,1,:,4) = ycomb3(:,:,1,:,4) + sum(sum(ycomb3(:,:,2:10,:,[10,15]),3),5);
            ycomb3_S(:,:,1,:,5) = ycomb3(:,:,1,:,5) + sum(sum(ycomb3(:,:,2:10,:,[11,16:20]),3),5);
            R0_testing(i,j) = R0(70,real(reshape(ycomb3_S(find(TT2_treat3>70,1),:,:,:,1:20),3*10*9*20,1)));
            for year_elim = 1:Run-1
                incred_test(i,year_elim) = (inc2015 - sum(sum(sum(ycomb3_noage(find(TT2_treat3>=Tin+year_elim,1)+1:find(TT2_treat3>=Tin+year_elim+1,1),:,:,27)))))/...
                    (TT2_treat3(find(TT2_treat3>=Tin+year_elim+1,1))-TT2_treat3(find(TT2_treat3>=Tin+year_elim,1))) / inc2015;
            end
            if incred_test(i,end) >= 0.8
                year_elim_test(i,j) = 1950+Tin+find(incred_test(i,:)>=0.8,1);
            else
                year_elim_test(i,j) = -1;
            end
            if R0_testing(i,j) < 1 && j > 1
                if R0_testing(i,1) > 1
                    infect_R0(i) = interp1(R0_testing(i,1:j), infect_red(1:j), 1);
                else
                    infect_R0(i) = -1;
                end
                break;                
            end
        end
    end
    R0_scens(s,4:7) = R0_testing';
    
    [cost_treat3,qaly_treat3,life_exp_treat3,tr3]=Costs_age(y2_treat3_noage(1:t_val_treat3(2)-t_val_treat3(1)+1,:,:,:),TT_3);
    tr3_(1) = sum(sum(sum(y2_treat3_noage(1:t_val_treat3(3)-t_val_treat3(1)+1,1,:,[23,25])))); tr3_(2) = sum(sum(sum(y2_treat3_noage(1:t_val_treat3(3)-t_val_treat3(1)+1,1,:,[23,25])))); tr3_(3) = sum(sum(sum(y2_treat3_noage(1:t_val_treat3(3)-t_val_treat3(1)+1,3,:,[23,25])))); tr3_(4) = sum(sum(sum(sum(y2_treat3_noage(1:t_val_treat3(3)-t_val_treat3(1)+1,:,:,[23,25])))));
    summary(3,:,s)=[cost_treat3,qaly_treat3,tr3(4),tr3_(4),...
        sum(sum(sum(ycomb3_noage(find(TT2_treat3>=80,1)+1:find(TT2_treat3>=81,1),:,:,27))))/(TT2_treat3(find(TT2_treat3>=81,1))-TT2_treat3(find(TT2_treat3>=80,1))),... %Incidence in 2030
        (sum(sum(ycomb3_noage(find(TT2_treat3>=81,1),1:3,:,22)))-sum(sum(ycomb3_noage(find(TT2_treat3>=80,1),1:3,:,22))))...
        /(TT2_treat3(find(TT2_treat3>=81,1))-TT2_treat3(find(TT2_treat3>=80,1))),... %Liver related deaths in 2030
        100*sum(sum(ycomb3_noage(t_val_treat3(1)-1,1,:,[6,12:20])))./sum(sum(ycomb3_noage(t_val_treat3(1),1,:,1:20))),... %Prevalence among PWID in 2015
        100*sum(sum(ycomb3_noage(t_val_treat3(3),1,:,[6,12:20])))./sum(sum(ycomb3_noage(t_val_treat3(3),1,:,1:20))),...%Prevalence among PWID in 2030
        sum(sum(ycomb3_noage(t_val_treat3(3),1:3,:,22))),... %Liver related deaths 2015 to 2030
        100*sum(sum(ycomb3_noage(find(TT2_treat3>=69,1),1,:,[6,12:20])))./sum(sum(ycomb3_noage(find(TT2_treat3>=69,1),1,:,1:20))),... %Prevalence among PWID in 2019
        100*sum(sum(ycomb3_noage(find(TT2_treat3>=70,1),1,:,[6,12:20])))./sum(sum(ycomb3_noage(find(TT2_treat3>=70,1),1,:,1:20))),... %Prevalence among PWID in 2020
        100*sum(sum(ycomb3_noage(find(TT2_treat3>=75,1),1,:,[6,12:20])))./sum(sum(ycomb3_noage(find(TT2_treat3>=75,1),1,:,1:20)))]; %Prevalence among PWID in 2025
    
    
    %%  Scenario 5: All for WHO
    scenario = 'current'
    %p_complete=1;
    %clear deathred death_scaleup incred testing_scaleup death_scale
    if strcmp(scenario,'WHO') == 1
        range = [1/c1_chronicPWID_base]; % Estimate of testing required
        range_death = fliplr(0.1:0.5:2.2);
        range_treat = 10:20:100;
        target_late = 0;
        %inc2015 = sum(sum(sum(ycomb_noage(find(TT2>=65,1)+1:find(TT2>=66,1),:,:,27))))/(TT2(find(TT2>=66,1))-TT2(find(TT2>=65,1)));
        %death2015 = (sum(sum(ycomb_noage(find(TT2>=66,1),1:3,:,22)))-sum(sum(ycomb_noage(find(TT2>=65,1),1:3,:,22))))/(TT2(find(TT2>=66,1))-TT2(find(TT2>=65,1)));
        for i =1:length(range)
            i
            c1_chronicPWID = range(i) * c1_chronicPWID_base;
            if c2_PWID < 100 c2_PWID = range(i) * c2_PWID_base; end
            for j = 1:length(range_death)
                j
                r_S4death = 1/(-1/log(1-(0.138+0.605)/2)-1/log(1-0.02 * range_death(j)));
                [TT2_treat4,y2_treat4]=DE_track_age(Run,y1_end,TT1,treat);
                ycomb4=[y1(1:end-1,:,:,:,:);y2_treat4];
                y2_treat4_noage=reshape(sum(y2_treat4(:,:,:,:,:),4),length(y2_treat4(:,1,1,1,1)),3,10,27+6);
                ycomb4_noage=reshape(sum(ycomb4(:,:,:,:,:),4),length(ycomb4(:,1,1,1,1)),3,10,27+6);
                deathred(i,j) = (death2015 - (sum(sum(ycomb4_noage(find(TT2_treat4>=81,1),1:3,:,22)))-sum(sum(ycomb4_noage(find(TT2_treat4>=80,1),1:3,:,22))))/(TT2_treat4(find(TT2_treat4>=81,1))-TT2_treat4(find(TT2_treat4>=80,1)))) / death2015;
                if deathred(i,j) > target_death && j > 1
                    break;
                end
            end
            if deathred(i,1) > target_death
                death_scaleup(i) = 1;
            else
                death_scaleup(i) = interp1(deathred(i,1:j), range_death(1:j), target_death);
            end
            r_S4death = 1/(-1/log(1-(0.138+0.605)/2)-1/log(1-0.02 * death_scaleup(i)));
            [TT2_treat4,y2_treat4]=DE_track_age(Run,y1_end,TT1,treat);
            ycomb4=[y1(1:end-1,:,:,:,:);y2_treat4];
            y2_treat4_noage=reshape(sum(y2_treat4(:,:,:,:,:),4),length(y2_treat4(:,1,1,1,1)),3,10,27+6);
            ycomb4_noage=reshape(sum(ycomb4(:,:,:,:,:),4),length(ycomb4(:,1,1,1,1)),3,10,27+6);
            incred(i) = (inc2015 - sum(sum(sum(ycomb4_noage(find(TT2_treat4>=80,1)+1:find(TT2_treat4>=81,1),:,:,27))))/(TT2_treat4(find(TT2_treat4>=81,1))-TT2_treat4(find(TT2_treat4>=80,1)))) / inc2015;
            if incred(i) > target_inc && i > 1
                break;
            end
        end
        
        testing_scaleup = interp1(incred, range(1:i), target_inc);
        death_scale = death_scaleup(find(range>=testing_scaleup,1));
        r_S4death = 1/(-1/log(1-(0.138+0.605)/2)-1/log(1-0.02 * death_scale));
        c1_chronicPWID = testing_scaleup * c1_chronicPWID_base;
        if c2_PWID < 100 c2_PWID = testing_scaleup * c2_PWID_base; end
    end
    
    %clear deathred death_scaleup incred testing_scaleup death_scale
    if strcmp(scenario,'WHO1') == 1
        range_death = fliplr(0.1:0.5:2.2);
        for j =1:length(range_death)
            %temp_rate = (c1_chronicPWID_base - 0.31*1) / (1-0.31);
            %c1_chronicPWID = 0.67*1 + (1-0.67)*temp_rate;
            %if c2_PWID < 100 c2_PWID = range(i) * c2_PWID_base; end
            r_S4death = 1/(-1/log(1-(0.138+0.605)/2)-1/log(1-0.02 * range_death(j)));
            [TT2_treat4,y2_treat4]=DE_track_age(Run,y1_end,TT1,treat);
            ycomb4=[y1(1:end-1,:,:,:,:);y2_treat4];
            y2_treat4_noage=reshape(sum(y2_treat4(:,:,:,:,:),4),length(y2_treat4(:,1,1,1,1)),3,10,27+6);
            ycomb4_noage=reshape(sum(ycomb4(:,:,:,:,:),4),length(ycomb4(:,1,1,1,1)),3,10,27+6);
            deathred(j) = (death2015 - (sum(sum(ycomb4_noage(find(TT2_treat4>=81,1),1:3,:,22)))-sum(sum(ycomb4_noage(find(TT2_treat4>=80,1),1:3,:,22))))/(TT2_treat4(find(TT2_treat4>=81,1))-TT2_treat4(find(TT2_treat4>=80,1)))) / death2015;
            if deathred(j) > target_death && j > 1
                break;
            end
        end
        if deathred(1) > target_death
            death_scale = 1;
        else
            death_scale = interp1(deathred(1:j), range_death(1:j), target_death);
        end
        r_S4death = 1/(-1/log(1-(0.138+0.605)/2)-1/log(1-0.02 * death_scale));
        [TT2_treat4,y2_treat4]=DE_track_age(Run,y1_end,TT1,treat);
        ycomb4=[y1(1:end-1,:,:,:,:);y2_treat4];
        y2_treat4_noage=reshape(sum(y2_treat4(:,:,:,:,:),4),length(y2_treat4(:,1,1,1,1)),3,10,27+6);
        ycomb4_noage=reshape(sum(ycomb4(:,:,:,:,:),4),length(ycomb4(:,1,1,1,1)),3,10,27+6);
    end
    
    
    [TT2_treat4,y2_treat4]=DE_track_age(Run,y1_end,TT1,treat);
    %p_complete=0.9; r_S4death = 1/(-1/log(1-(0.138+0.605)/2)-1/log(1-0.01));
    
    ycomb4=[y1(1:end-1,:,:,:,:);y2_treat4];
    y2_treat4_noage=reshape(sum(y2_treat4(:,:,:,:,:),4),length(y2_treat4(:,1,1,1,1)),3,10,27+6);
    ycomb4_noage=reshape(sum(ycomb4(:,:,:,:,:),4),length(ycomb4(:,1,1,1,1)),3,10,27+6);
    
    Tint=Tin+[0,3,14]; t_val_treat4=zeros(1,length(Tint));
    for l=1:length(Tint) %Find times corresponding to years since intervention
        t_val_treat4(l)=find(TT2_treat4>=(Tint(l)),1);
    end
    TT_4=TT2_treat4(t_val_treat4(1):t_val_treat4(2))-Tin;
    
    
    [cost_treat4,qaly_treat4,life_exp_treat4,tr4]=Costs_age(y2_treat4_noage(1:t_val_treat4(2)-t_val_treat4(1)+1,:,:,:),TT_4);
    tr4_(1) = sum(sum(sum(y2_treat4_noage(1:t_val_treat4(3)-t_val_treat4(1)+1,1,:,[23,25])))); tr4_(2) = sum(sum(sum(y2_treat4_noage(1:t_val_treat4(3)-t_val_treat4(1)+1,1,:,[23,25])))); tr4_(3) = sum(sum(sum(y2_treat4_noage(1:t_val_treat4(3)-t_val_treat4(1)+1,3,:,[23,25])))); tr4_(4) = sum(sum(sum(sum(y2_treat4_noage(1:t_val_treat4(3)-t_val_treat4(1)+1,:,:,[23,25])))));
    summary(4,:,s)=[cost_treat4,qaly_treat4,tr4(4),tr4_(4),...
        sum(sum(sum(ycomb4_noage(find(TT2_treat4>=80,1)+1:find(TT2_treat4>=81,1),:,:,27))))/(TT2_treat4(find(TT2_treat4>=81,1))-TT2_treat4(find(TT2_treat4>=80,1))),... %Incidence in 2030
        (sum(sum(ycomb4_noage(find(TT2_treat4>=81,1),1:3,:,22)))-sum(sum(ycomb4_noage(find(TT2_treat4>=80,1),1:3,:,22))))...
        /(TT2_treat4(find(TT2_treat4>=81,1))-TT2_treat4(find(TT2_treat4>=80,1))),... %Liver related deaths in 2030
        100*sum(sum(ycomb4_noage(t_val_treat4(1)-1,1,:,[6,12:20])))./sum(sum(ycomb4_noage(t_val_treat4(1),1,:,1:20))),... %Prevalence among PWID in 2015
        100*sum(sum(ycomb4_noage(t_val_treat4(3),1,:,[6,12:20])))./sum(sum(ycomb4_noage(t_val_treat4(3),1,:,1:20))),...%Prevalence among PWID in 2030
        sum(sum(ycomb4_noage(t_val_treat4(3),1:3,:,22))),... %Liver related deaths 2015 to 2030
        100*sum(sum(ycomb4_noage(find(TT2_treat4>=69,1),1,:,[6,12:20])))./sum(sum(ycomb4_noage(find(TT2_treat4>=69,1),1,:,1:20))),... %Prevalence among PWID in 2019
        100*sum(sum(ycomb4_noage(find(TT2_treat4>=70,1),1,:,[6,12:20])))./sum(sum(ycomb4_noage(find(TT2_treat4>=70,1),1,:,1:20))),... %Prevalence among PWID in 2020
        100*sum(sum(ycomb4_noage(find(TT2_treat4>=75,1),1,:,[6,12:20])))./sum(sum(ycomb4_noage(find(TT2_treat4>=75,1),1,:,1:20)))]; %Prevalence among PWID in 2025
    
    %% Collect outputs
    
    years=60:1:81;
    for t=2:length(years)
        inc_year_sens(:,t-1,s)=[sum(sum(sum(ycomb_noage(find(TT2>=years(t-1),1)+1:find(TT2>=years(t),1),:,:,27))));...
            sum(sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1)+1:find(TT2_treat>=years(t),1),:,:,27))));...
            sum(sum(sum(ycomb3_noage(find(TT2_treat3>=years(t-1),1)+1:find(TT2_treat3>=years(t),1),:,:,27))));...
            sum(sum(sum(ycomb4_noage(find(TT2_treat4>=years(t-1),1)+1:find(TT2_treat4>=years(t),1),:,:,27))));...
            sum(sum(sum(ycomb5_noage(find(TT2_treat5>=years(t-1),1)+1:find(TT2_treat5>=years(t),1),:,:,27))));...
            sum(sum(sum(ycomb6_noage(find(TT2_treat6>=years(t-1),1)+1:find(TT2_treat6>=years(t),1),:,:,27))))]...
            ./[TT2(find(TT2>=years(t),1))-TT2(find(TT2>=years(t-1),1)+1);TT2_treat(find(TT2_treat>=years(t),1))-TT2_treat(find(TT2_treat>=years(t-1),1)+1);...
            TT2_treat3(find(TT2_treat3>=years(t),1))-TT2_treat3(find(TT2_treat3>=years(t-1),1)+1);TT2_treat4(find(TT2_treat4>=years(t),1))-TT2_treat4(find(TT2_treat4>=years(t-1),1)+1);...
            TT2_treat5(find(TT2_treat5>=years(t),1))-TT2_treat5(find(TT2_treat5>=years(t-1),1)+1);TT2_treat6(find(TT2_treat6>=years(t),1))-TT2_treat6(find(TT2_treat6>=years(t-1),1)+1)];
        
        death_year_sens(:,t-1,s)=[sum(sum(ycomb_noage(find(TT2>=years(t),1),:,:,22))) - sum(sum(ycomb_noage(find(TT2>=years(t-1),1),:,:,22)));...
            sum(sum(ycomb2_noage(find(TT2_treat>=years(t),1),:,:,22))) - sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1),:,:,22)));...
            sum(sum(ycomb3_noage(find(TT2_treat3>=years(t),1),:,:,22))) - sum(sum(ycomb3_noage(find(TT2_treat3>=years(t-1),1),:,:,22)));...
            sum(sum(ycomb4_noage(find(TT2_treat4>=years(t),1),:,:,22))) - sum(sum(ycomb4_noage(find(TT2_treat4>=years(t-1),1),:,:,22)));...
            sum(sum(ycomb5_noage(find(TT2_treat5>=years(t),1),:,:,22))) - sum(sum(ycomb5_noage(find(TT2_treat5>=years(t-1),1),:,:,22)));...
            sum(sum(ycomb6_noage(find(TT2_treat6>=years(t),1),:,:,22))) - sum(sum(ycomb6_noage(find(TT2_treat6>=years(t-1),1),:,:,22)))]...
            ./...
            [TT2(find(TT2>=years(t),1))-TT2(find(TT2>=years(t-1),1));...
            TT2_treat(find(TT2_treat>=years(t),1))-TT2_treat(find(TT2_treat>=years(t-1),1));...
            TT2_treat3(find(TT2_treat3>=years(t),1))-TT2_treat3(find(TT2_treat3>=years(t-1),1));...
            TT2_treat4(find(TT2_treat4>=years(t),1))-TT2_treat4(find(TT2_treat4>=years(t-1),1));...
            TT2_treat5(find(TT2_treat5>=years(t),1))-TT2_treat5(find(TT2_treat5>=years(t-1),1));...
            TT2_treat6(find(TT2_treat6>=years(t),1))-TT2_treat6(find(TT2_treat6>=years(t-1),1))];
        
        treatment_alloc_sens(:,t-1,s)=[sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1):find(TT2_treat>years(t),1)-1,1,:,23)));...% PWID early liver disease treatments
            sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1):find(TT2_treat>years(t),1)-1,1,:,25)));... PWID late liver disease treatments
            sum(sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1):find(TT2_treat>years(t),1)-1,2:3,:,23))));... %Former PWID early liver disease stage
            sum(sum(sum(ycomb2_noage(find(TT2_treat>=years(t-1),1):find(TT2_treat>years(t),1)-1,2:3,:,25))))];%... %Former PWID late liver disease stage
        %./repmat(TT2_treatlate(find(TT2_treatlate>=years(t),1)-1)-TT2_treatlate(find(TT2_treatlate>=years(t-1),1)),4,1);
    end
    years2=50:1:81;
    for t=2:length(years2)
        HCC_year_sens(:,t-1,s)=[sum(sum(sum(ycomb_noage(find(TT2>=years2(t-1),1)+1:find(TT2>=years2(t),1),:,:,24))));...
            sum(sum(sum(ycomb2_noage(find(TT2_treat>=years2(t-1),1)+1:find(TT2_treat>=years2(t),1),:,:,24))));...
            sum(sum(sum(ycomb3_noage(find(TT2_treat3>=years2(t-1),1)+1:find(TT2_treat3>=years2(t),1),:,:,24))));...
            sum(sum(sum(ycomb4_noage(find(TT2_treat4>=years2(t-1),1)+1:find(TT2_treat4>=years2(t),1),:,:,24))));...
            sum(sum(sum(ycomb5_noage(find(TT2_treat5>=years2(t-1),1)+1:find(TT2_treat5>=years2(t),1),:,:,24))));...
            sum(sum(sum(ycomb6_noage(find(TT2_treat6>=years2(t-1),1)+1:find(TT2_treat6>=years2(t),1),:,:,24))))]...
            ./[TT2(find(TT2>=years2(t),1))-TT2(find(TT2>=years2(t-1),1)+1);TT2_treat(find(TT2_treat>=years2(t),1))-TT2_treat(find(TT2_treat>=years2(t-1),1)+1);...
            TT2_treat3(find(TT2_treat3>=years2(t),1))-TT2_treat3(find(TT2_treat3>=years2(t-1),1)+1);TT2_treat4(find(TT2_treat4>=years2(t),1))-TT2_treat4(find(TT2_treat4>=years2(t-1),1)+1);...
            TT2_treat5(find(TT2_treat5>=years2(t),1))-TT2_treat5(find(TT2_treat5>=years2(t-1),1)+1);TT2_treat6(find(TT2_treat6>=years2(t),1))-TT2_treat6(find(TT2_treat6>=years2(t-1),1)+1)];
        
        diagnosed_year_sens(:,t-1,s)=[sum(sum(sum(ycomb_noage(find(TT2>=years2(t-1),1)+1:find(TT2>=years2(t),1),:,:,29))));...
            sum(sum(sum(ycomb2_noage(find(TT2_treat>=years2(t-1),1)+1:find(TT2_treat>=years2(t),1),:,:,29))));...
            sum(sum(sum(ycomb3_noage(find(TT2_treat3>=years2(t-1),1)+1:find(TT2_treat3>=years2(t),1),:,:,29))));...
            sum(sum(sum(ycomb4_noage(find(TT2_treat4>=years2(t-1),1)+1:find(TT2_treat4>=years2(t),1),:,:,29))));...
            sum(sum(sum(ycomb5_noage(find(TT2_treat5>=years2(t-1),1)+1:find(TT2_treat5>=years2(t),1),:,:,29))));...
            sum(sum(sum(ycomb6_noage(find(TT2_treat6>=years2(t-1),1)+1:find(TT2_treat6>=years2(t),1),:,:,29))))]...
            ./[TT2(find(TT2>=years2(t),1))-TT2(find(TT2>=years2(t-1),1)+1);TT2_treat(find(TT2_treat>=years2(t),1))-TT2_treat(find(TT2_treat>=years2(t-1),1)+1);...
            TT2_treat3(find(TT2_treat3>=years2(t),1))-TT2_treat3(find(TT2_treat3>=years2(t-1),1)+1);TT2_treat4(find(TT2_treat4>=years2(t),1))-TT2_treat4(find(TT2_treat4>=years2(t-1),1)+1);...
            TT2_treat5(find(TT2_treat5>=years2(t),1))-TT2_treat5(find(TT2_treat5>=years2(t-1),1)+1);TT2_treat6(find(TT2_treat6>=years2(t),1))-TT2_treat6(find(TT2_treat6>=years2(t-1),1)+1)];
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

%Treatment numbers
tot_infected_PWID=sum(sum(ycomb_noage(t_val(1),1,:,[6,12:20])));
tot_PWID_adv=sum(sum(ycomb_noage(t_val(1),1,:,[15:20])));
tot_former_adv=sum(sum(sum(ycomb_noage(t_val(1),2:3,:,[15:20]))));
if target_late<0
    target_late=sum(sum(sum(ycomb_noage(t_val(1),:,:,15:20))))/sum(sum(sum(ycomb_noage(t_val(1),:,:,12:20))));
else
    intersection=tot_PWID_adv/(target_late*(tot_PWID_adv + tot_former_adv) + (1-target_late)*tot_infected_PWID); %The proportion of treatments to be allocated to infected PWID with late liver disease
end


inc_year=mean(inc_year_sens(:,:,:)./repmat(inc_year_sens(1,1,:),6,length(years)-1),3);
inc_year(5,:)=mean(inc_year_sens(5,:,treat_scaleup(1,:)>=0)./repmat(inc_year_sens(1,1,treat_scaleup(1,:)>=0),1,length(years)-1),3);

%inc_year = mean(inc_year_sens(:,:,:),3);
inc_year_std=prctile((inc_year_sens(:,:,:)./repmat(inc_year_sens(1,1,:),6,length(years)-1)),[5,95],3);
inc_year_std(5,:,:)=prctile((inc_year_sens(5,:,treat_scaleup(1,:)>=0)./repmat(inc_year_sens(1,1,treat_scaleup(1,:)>=0),1,length(years)-1)),[5,95],3);
deaths15_30=mean(sum(death_year_sens(:,[1:15],:),2),3);
deaths15_30_std=prctile(sum(death_year_sens(:,[1:15],:),2),[5,95],3);
death_year=mean(death_year_sens(:,:,:),3);
death_year_std=prctile((death_year_sens(:,:,:)./repmat(death_year_sens(1,1,:),6,length(years)-1)),[5,95],3);
treatment_alloc=mean(treatment_alloc_sens,3);
treatment_alloc_std=std(treatment_alloc_sens,0,3);
HCC_year=mean(HCC_year_sens(:,:,:),3);
diagnosed_year=mean(diagnosed_year_sens(:,:,:),3);
treat_scaleup_summary(:,1) = mean(treat_scaleup,2); treat_scaleup_summary(1,1) = mean(treat_scaleup(1,treat_scaleup(1,:)>=0),2);
treat_scaleup_summary(:,2:3) = prctile(treat_scaleup,[5,96],2); treat_scaleup_summary(1,2:3) = prctile(treat_scaleup(1,treat_scaleup(1,:)>=0),[5,95],2);
R0_summary(1,:) = mean(R0_scens,1); R0_summary(2:3,:) = prctile(R0_scens, [5,95],1);


for i = 1:length(Total_HCV_sens(:,1,1))
    for j = 1:6
        if i==5
            Total_HCV(i,j,:) = prctile(Total_HCV_sens(i,j,treat_scaleup(1,:)>=0),[5,95]);
        else
            Total_HCV(i,j,:) = prctile(Total_HCV_sens(i,j,:),[5,95]);
        end
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
margins(:,:,1)=mean(summary(:,:,:),3); margins(5,:,1)=mean(summary(5,:,treat_scaleup(1,:)>=0),3);
margins2=zeros(6,12,3); 
margins2(:,:,1)=mean(summary2(:,:,:),3); margins2(5,:,1)=mean(summary2(5,:,treat_scaleup(1,:)>=0),3);
total_treat_summary = zeros(6,4,3); 
total_treat_summary(:,:,1) = mean(total_treat(:,:,:),3); total_treat_summary(5,:,1) = mean(total_treat(5,:,treat_scaleup(1,:)>=0),3);
total_treat_2030_summary = zeros(6,4,3); 
total_treat_2030_summary(:,:,1) = mean(total_treat_2030(:,:,:),3); total_treat_2030_summary(5,:,1) = mean(total_treat_2030(5,:,treat_scaleup(1,:)>=0),3);

for i=1:6
    for j=1:12
        if i==5
            margins(i,j,2)=prctile(summary(i,j,treat_scaleup(1,:)>=0),5); %Lower bound
            margins(i,j,3)=prctile(summary(i,j,treat_scaleup(1,:)>=0),95); %Upper bound
        
            margins2(i,j,2)=prctile(summary2(i,j,treat_scaleup(1,:)>=0),5); %Lower bound
            margins2(i,j,3)=prctile(summary2(i,j,treat_scaleup(1,:)>=0),95); %Upper bound
        else
            margins(i,j,2)=prctile(summary(i,j,:),5); %Lower bound
            margins(i,j,3)=prctile(summary(i,j,:),95); %Upper bound
        
            margins2(i,j,2)=prctile(summary2(i,j,:),5); %Lower bound
            margins2(i,j,3)=prctile(summary2(i,j,:),95); %Upper bound
        end
    end
    for j=1:4
        if i==5
            total_treat_summary(i,j,2)=prctile(total_treat(i,j,treat_scaleup(1,:)>=0),5); %Lower bound
            total_treat_summary(i,j,3)=prctile(total_treat(i,j,treat_scaleup(1,:)>=0),95); %Upper bound
            total_treat_2030_summary(i,j,2)=prctile(total_treat_2030(i,j,treat_scaleup(1,:)>=0),5); %Lower bound
            total_treat_2030_summary(i,j,3)=prctile(total_treat_2030(i,j,treat_scaleup(1,:)>=0),95); %Upper bound
        else
            total_treat_summary(i,j,2)=prctile(total_treat(i,j,:),5); %Lower bound
            total_treat_summary(i,j,3)=prctile(total_treat(i,j,:),95); %Upper bound
            total_treat_2030_summary(i,j,2)=prctile(total_treat_2030(i,j,:),5); %Lower bound
            total_treat_2030_summary(i,j,3)=prctile(total_treat_2030(i,j,:),95); %Upper bound
        end
    end
end

testing = 1/c1_chronicPWID;
save(filename, ...
    'summary', 'summary2', 'margins', 'margins2', 'inc_year', 'death_year',...
    'testing', 'total_treat_summary','total_treat_2030_summary',...
    'Total_HCV_sens', 'treat_scaleup', 'year_elim_test','R0_scens', 'R0_testing')


alpha = alpha_old;
infect=infect_base;
c1_chronicPWID = c1_chronicPWID_base;
c2_PWID = c2_PWID_base;
c3_PWID = c3_PWID_base;
c4_PWID = c4_PWID_base;
c1_chronicformerPWID = c1_chronicformerPWID_base;
c2_formerPWID = c2_formerPWID_base;
c3_formerPWID = c3_formerPWID_base;
c4_formerPWID = c4_formerPWID_base;
c1_chronic = c1_chronic_base;
c2 = c2_base;
c3 = c3_base;
c4 = c4_base;
imp1 = imp1_base;
imp2 = imp2_base;
imp3 = imp3_base;
imp4 = imp4_base;
imp5 = imp5_base;
imp6 = imp6_base;
imp7 = imp7_base;
imp8 = imp8_base;
imp9 = imp9_base;
r_S4death = r_S4death_base;


end

