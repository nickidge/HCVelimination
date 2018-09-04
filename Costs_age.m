function[c,DALY,life,tr,c_HR]=Costs_age(y,TT,y_large)
global Q_sens q_svr q_svr_PWID q_treat sens c_daa discount care cascade_scale_time scenario prop_test progression ...
    num_pops num_cascade num_age num_intervention 

%% Costs
[costing0,costing0_label,costing0_raw] = xlsread('Template.xlsx','costs'); 
for i = 1:length(costing0_raw(2:end,1))
    str = costing0_raw(i+1,1);
    p = strcmpi(str,costing0_raw(:,1));
    evalc(sprintf('%s=%s',cell2mat(str),num2str(cell2mat(costing0_raw(p,2)))));
end

% disease stage costs
c_svr0=[c_svr0__a, c_svr0__b];
c_svr_pa=[c_svr_pa__a,c_svr_pa__b,c_svr_pa__c]; %For successful treatment at F0-F2, F3 and F4+ stages
c_cascade3=[c_cascade3__a,c_cascade3__b];
c_cascade4=[c_cascade4__a,c_cascade4__b,c_cascade4__c];
c_daa = [c_daa__a,c_daa__b];

%% QALYS
if sens==1
    [health0,health0_label,health0_raw] = xlsread('Template.xlsx','health_utilities');
    for i = 1:length(health0_raw(2:end,1))
        str = health0_raw(i+1,1);
        p = strcmpi(str,health0_raw(:,1));
        evalc(sprintf('%s=%s',cell2mat(str),num2str(cell2mat(health0_raw(p,2)))));
    end
    
    q_svr_PWID=[q_svr_PWID__a,q_svr_PWID__b,q_svr_PWID__c]; %After successful treatment from F012, or treatment from F3, or treatment from F4 onwards
    q_svr=[q_svr__a,q_svr__b,q_svr__c];
    q_treat=([q_treat__a,q_treat__b,q_treat__c]+q_svr)/2; %During treatment from F012; F3; F4 or later stages
end

%% QALYS: Sensitivity analysis
if sens>1
    q_S_PWID=Q_sens(1);
    q_S=Q_sens(2);
    q_A=Q_sens(3);
    q_F012=Q_sens(4);
    q_F3=Q_sens(5);
    q_F4=Q_sens(6);
    q_DC=Q_sens(7);
    q_HCC=Q_sens(8);
    q_LT1=Q_sens(9);
    q_LT2=Q_sens(10);
end

%% Costs / QALY among PWID

C1a=trapz(TT,sum(sum(y(:,1,:,1:3),3),4).*exp(-discount.*TT))*c_svr_pa(1);
C1b=trapz(TT,sum(y(:,1,:,4),3).*exp(-discount.*TT))*c_svr_pa(2);
C1c=trapz(TT,sum(y(:,1,:,5),3).*exp(-discount.*TT))*c_svr_pa(3);
C2a=sum((sum(y(:,1,:,23),3)*(c_daa(1)+c_PWID)+sum(y(:,1,:,23),3)*0.9*c_svr0(1)).*exp(-discount.*TT)); %Cost of treatment from the F0-F3 stages
C2b=sum((sum(y(:,1,:,25),3)*(c_daa(2)+c_PWID)+sum(y(:,1,:,25),3)*0.9*c_svr0(2)).*exp(-discount.*TT)); %Cost of treatment from the F4 stage onwards
C3=trapz(TT,sum(y(:,1,1:10,6),3).*exp(-discount.*TT))*c_A; %Cost of diagnosis (measured in acute state). ONLY IF DIAGNOSED
C4=trapz(TT,sum(y(:,1,1:10,12),3).*exp(-discount.*TT))*c_F012; %Cost of F0 state (integral under curve. ONLY IF DIAGNOSED
C5=trapz(TT,sum(y(:,1,1:10,13),3).*exp(-discount.*TT))*c_F012;
C6=trapz(TT,sum(y(:,1,1:10,14),3).*exp(-discount.*TT))*c_F012;
C7=trapz(TT,sum(y(:,1,1:10,15),3).*exp(-discount.*TT))*c_F3;
C8=trapz(TT,sum(y(:,1,:,16),3).*exp(-discount.*TT))*c_F4;
C9=trapz(TT,sum(y(:,1,:,17),3).*exp(-discount.*TT))*c_DC;
C10=trapz(TT,sum(y(:,1,:,18),3).*exp(-discount.*TT))*c_HCC_pa;

C_F4transfer=sum(sum(y(:,1,:,21),3)*c_F4_0.*exp(-discount.*TT)); %One off cost of transition to F4 state
C_HCCtransfer=sum(sum(y(:,1,:,24),3)*c_HCC_0.*exp(-discount.*TT)); %One off cost of transition to HCC state
C_LT=sum(sum(y(:,1,:,26),3)*c_LT.*exp(-discount.*TT)); %One off cost of liver transplant


c_PWIDtot=(C1a+C1b+C1c+C2a+C2b+C3+C4+C5+C6+C7+C8+C9+C10+C_F4transfer+C_HCCtransfer+C_LT);

Q1a=trapz(TT,sum(sum(y(:,1,:,1:3),3),4).*exp(-discount.*TT))*q_svr_PWID(1);
Q1b=trapz(TT,sum(y(:,1,:,4),3).*exp(-discount.*TT))*q_svr_PWID(2);
Q1c=trapz(TT,sum(y(:,1,:,5),3).*exp(-discount.*TT))*q_svr_PWID(3);
Q3=trapz(TT,sum(y(:,1,:,6),3).*exp(-discount.*TT))*q_A; %QALYs in acute state (integral under curve)
Q4a=trapz(TT,sum(sum(y(:,1,:,7:9),3),4).*exp(-discount.*TT))*q_treat(1);
Q4b=trapz(TT,sum(y(:,1,:,10),3).*exp(-discount.*TT))*q_treat(2);
Q4c=trapz(TT,sum(y(:,1,:,11),3).*exp(-discount.*TT))*q_treat(3);
Q5=trapz(TT,sum(y(:,1,:,12),3).*exp(-discount.*TT))*q_F012; %QALYs in F0 state (integral under curve)
Q6=trapz(TT,sum(y(:,1,:,13),3).*exp(-discount.*TT))*q_F012;
Q7=trapz(TT,sum(y(:,1,:,14),3).*exp(-discount.*TT))*q_F012;
Q8=trapz(TT,sum(y(:,1,:,15),3).*exp(-discount.*TT))*q_F3;
Q9=trapz(TT,sum(y(:,1,:,16),3).*exp(-discount.*TT))*q_F4;
Q10=trapz(TT,sum(y(:,1,:,17),3).*exp(-discount.*TT))*q_DC;
Q11=trapz(TT,sum(y(:,1,:,18),3).*exp(-discount.*TT))*q_HCC;
Q12=trapz(TT,sum(y(:,1,:,19),3).*exp(-discount.*TT))*q_LT1;
Q13=trapz(TT,sum(y(:,1,:,20),3).*exp(-discount.*TT))*q_LT2;

for i=2:size(y,1) %Make sure these are the additional treatments/HCC transfers/LTs/liver transplants/incidence, rather than total
    y(i,:,:,22)=max(y(i,:,:,22)-y(i-1,:,:,22),0);
end
DALY1=trapz(TT,sum(y(:,1,:,17),3).*exp(-discount.*TT))*DALY_DC;
DALY2=trapz(TT,sum(y(:,1,:,18),3).*exp(-discount.*TT))*DALY_HCC;
DALY3=trapz(TT,sum(y(:,1,:,22),3).*exp(-discount.*TT))*DALY_death;

DALY_PWID = (DALY1+DALY2+DALY3);

q_PWID=(Q1a+Q1b+Q1c+Q3+Q4a+Q4b+Q4c+Q5+Q6+Q7+Q8+Q9+Q10+Q11+Q12+Q13); 


%% Costs / QALY among former PWID and non injectors

C1a=trapz(TT,sum(sum(sum(y(:,2:3,:,1:3),2),3),4).*exp(-discount.*TT))*c_svr_pa(1);
C1b=trapz(TT,sum(sum(y(:,2:3,:,4),2),3).*exp(-discount.*TT))*c_svr_pa(2);
C1c=trapz(TT,sum(sum(y(:,2:3,:,5),2),3).*exp(-discount.*TT))*c_svr_pa(3);
C2a=sum((sum(sum(y(:,2:3,:,23),2),3)*(c_daa(1)+c_formerPWID)+sum(sum(y(:,2:3,:,23),2),3)*0.9*c_svr0(1)).*exp(-discount.*TT)); %Cost of treatment from the F0-F3 stages
C2b=sum((sum(sum(y(:,2:3,:,25),2),3)*(c_daa(2)+c_formerPWID)+sum(sum(y(:,2:3,:,25),2),3)*0.9*c_svr0(2)).*exp(-discount.*TT)); %Cost of treatment from the F4 stage onwards
C3=trapz(TT,sum(sum(y(:,2:3,1:10,6),2),3).*exp(-discount.*TT))*c_A; %Cost of diagnosis (measured in acute state). ONLY IF DIAGNOSED
C4=trapz(TT,sum(sum(y(:,2:3,1:10,12),2),3).*exp(-discount.*TT))*c_F012; %Cost of F0 state (integral under curve)
C5=trapz(TT,sum(sum(y(:,2:3,1:10,13),2),3).*exp(-discount.*TT))*c_F012;
C6=trapz(TT,sum(sum(y(:,2:3,1:10,14),2),3).*exp(-discount.*TT))*c_F012;
C7=trapz(TT,sum(sum(y(:,2:3,1:10,15),2),3).*exp(-discount.*TT))*c_F3;
C8=trapz(TT,sum(sum(y(:,2:3,:,16),2),3).*exp(-discount.*TT))*c_F4;
C9=trapz(TT,sum(sum(y(:,2:3,:,17),2),3).*exp(-discount.*TT))*c_DC;
C10=trapz(TT,sum(sum(y(:,2:3,:,18),2),3).*exp(-discount.*TT))*c_HCC_pa;

C_F4transfer=sum(sum(sum(y(:,2:3,:,21),2),3)*c_F4_0.*exp(-discount.*TT)); %One off cost of transition to F4 state
C_HCCtransfer=sum(sum(sum(y(:,2:3,:,24),2),3)*c_HCC_0.*exp(-discount.*TT)); %One off cost of transition to HCC state
C_LT=sum(sum(sum(y(:,2:3,:,26),2),3)*c_LT.*exp(-discount.*TT)); %One off cost of liver transplant

c_formertot=(C1a+C1b+C1c+C2a+C2b+C3+C4+C5+C6+C7+C8+C9+C10+C_F4transfer+C_HCCtransfer+C_LT); 

Q1a=trapz(TT,sum(sum(sum(y(:,2:3,:,1:3),2),3),4).*exp(-discount.*TT))*q_svr(1);
Q1b=trapz(TT,sum(sum(y(:,2:3,:,4),2),3).*exp(-discount.*TT))*q_svr(2);
Q1c=trapz(TT,sum(sum(y(:,2:3,:,5),2),3).*exp(-discount.*TT))*q_svr(3);
Q3=trapz(TT,sum(sum(y(:,2:3,:,6),2),3).*exp(-discount.*TT))*q_A; %QALYs in acute state (integral under curve)
Q4a=trapz(TT,sum(sum(sum(y(:,2:3,:,7:9),2),3),4).*exp(-discount.*TT))*q_treat(1);
Q4b=trapz(TT,sum(sum(y(:,2:3,:,10),2),3).*exp(-discount.*TT))*q_treat(2);
Q4c=trapz(TT,sum(sum(y(:,2:3,:,11),2),3).*exp(-discount.*TT))*q_treat(3);
Q5=trapz(TT,sum(sum(y(:,2:3,:,12),2),3).*exp(-discount.*TT))*q_F012; %QALYs in F0 state (integral under curve)
Q6=trapz(TT,sum(sum(y(:,2:3,:,13),2),3).*exp(-discount.*TT))*q_F012;
Q7=trapz(TT,sum(sum(y(:,2:3,:,14),2),3).*exp(-discount.*TT))*q_F012;
Q8=trapz(TT,sum(sum(y(:,2:3,:,15),2),3).*exp(-discount.*TT))*q_F3;
Q9=trapz(TT,sum(sum(y(:,2:3,:,16),2),3).*exp(-discount.*TT))*q_F4;
Q10=trapz(TT,sum(sum(y(:,2:3,:,17),2),3).*exp(-discount.*TT))*q_DC;
Q11=trapz(TT,sum(sum(y(:,2:3,:,18),2),3).*exp(-discount.*TT))*q_HCC;
Q12=trapz(TT,sum(sum(y(:,2:3,:,19),2),3).*exp(-discount.*TT))*q_LT1;
Q13=trapz(TT,sum(sum(y(:,2:3,:,20),2),3).*exp(-discount.*TT))*q_LT2;


DALY1=trapz(TT,sum(sum(y(:,2:3,:,17),3),2).*exp(-discount.*TT))*DALY_DC;
DALY2=trapz(TT,sum(sum(y(:,2:3,:,18),3),2).*exp(-discount.*TT))*DALY_HCC;
DALY3=trapz(TT,sum(sum(y(:,2:3,:,22),3),2).*exp(-discount.*TT))*DALY_death;

DALY_other = (DALY1+DALY2+DALY3);

q_former=(Q1a+Q1b+Q1c+Q3+Q4a+Q4b+Q4c+Q5+Q6+Q7+Q8+Q9+Q10+Q11+Q12+Q13); 

%% Cascade costs
if strcmp(scenario,'rapidRNA') == 1 || strcmp(scenario,'WHO') == 1 || strcmp(scenario,'WHO1') == 1 || strcmp(scenario,'current') == 1
    c_cascade2 = 163.90; % cost of rapid RNA test
    Cas_cost1=sum(sum(y(:,1,1,1:20),4)* prop_test * progression(1,1,2,1) * c_cascade2_ab.*exp(-discount.*TT)); % ab testing, for all PWID
    Cas_cost2=sum(sum(sum(y(:,:,:,29),3),2)*c_cascade2_RNA.*exp(-discount.*TT)); % i.e. RNA tests
    
elseif strcmp(scenario,'serum_HCVcAg') == 1
    Cas_cost1=sum(sum(y(:,1,1,1:20),4)* prop_test * progression(1,1,2,1) * c_cascade2_serum.*exp(-discount.*TT)); % i.e. serum tests only, for all PWID
    Cas_cost2 = 0;
elseif strcmp(scenario,'DBS_HCVcAg') == 1
    Cas_cost1=sum(sum(y(:,1,1,1:20),4)* prop_test * progression(1,1,2,1) * c_cascade2_DBS.*exp(-discount.*TT)); % i.e. DBS tests only, for all PWID
    Cas_cost2 = 0;
else
    Cas_cost1=sum(sum(y(:,1,1,1:20),4)* prop_test * progression(1,1,2,1) * c_cascade1.*exp(-discount.*TT)); % i.e. antibody tests
    Cas_cost2=sum(sum(sum(y(:,:,:,29),3),2)*c_cascade2.*exp(-discount.*TT)); % i.e RNA tests
end
Cas_cost3=sum(sum(sum(y(:,:,:,30),3),2)*c_cascade3(2).*exp(-discount.*TT)); % i.e assumes no genotype tests, just other workup

%Initial treatment
Cas_cost4=sum(sum(sum(y(:,:,6,32),3),2).*c_cascade4(3).*exp(-discount.*TT)); % late stage liver assessment required by specialist
if cascade_scale_time == 0 % scenarios where we only use the final GP:specialist ratio
    Cas_cost4=sum(sum(y(:,:,6,31),2).*(care(2) * c_cascade4(1) + (1-care(2)) * c_cascade4(3)).*exp(-discount.*TT)); % early stage first treatment
else
    t = min(1,TT / cascade_scale_time);
    Cas_cost4=sum(sum(y(:,:,6,31),2).* (...
        (1-t)*(care(1) * c_cascade4(1) + (1-care(1)) * c_cascade4(3)) + t * (care(2) * c_cascade4(1) + (1-care(2)) * c_cascade4(3))...
        ).*exp(-discount.*TT));
end
% retreatment
Cas_retreat=sum((sum(sum(y(:,:,8,[23,25]),4),2)*c_cascade4(3).*exp(-discount.*TT))); %Cost of workup for re-treatment from all liver disease stages (requiring specialist)

c_cascade=(Cas_cost1+Cas_cost2+Cas_cost3+Cas_cost4+Cas_retreat); 

%% Harm reduction
y2=squeeze(sum(sum(sum(y_large,4),6),7)); %Reshape to sum over age stratification, interventions, engagement and region
OST = trapz(TT,sum(sum(sum(y2(:,1,:,3:4,1:20),3),4),5).*exp(-discount.*TT))*c_OST;
NSP = trapz(TT,sum(sum(sum(y2(:,1,:,[2,4],1:20),3),4),5).*exp(-discount.*TT))*c_NSP;

c_HR = OST + NSP;

%% Total
c=c_PWIDtot+c_formertot+c_cascade+c_HR;
q=q_PWID+q_former;
life=trapz(TT,sum(sum(sum(y(:,:,:,1:20),4),3),2)) / sum(sum(sum(y(1,:,:,1:20))));
DALY = DALY_PWID+DALY_other;
%%treatments

tr(1) = sum(sum(sum(y(:,1,:,[23,25]))));
tr(2) = sum(sum(sum(y(:,2,:,[23,25]))));
tr(3) = sum(sum(sum(y(:,3,:,[23,25]))));
tr(4) = sum(sum(sum(sum(y(:,:,:,[23,25])))));

end

