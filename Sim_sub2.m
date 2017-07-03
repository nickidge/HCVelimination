function[ydot]=Sim_sub2(TT,y)
global mu_PWID mu_former exit_IDU r_relapse delta alpha p_complete omega infect total_PWID PWID0...
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2...
    c1_chronicPWID c2_PWID c3_PWID c4_PWID c5_PWID c6_PWID ...
    c1_chronicformerPWID c2_formerPWID c3_formerPWID c4_formerPWID c5_formerPWID c6_formerPWID ...
    c1_chronic c2 c3 c4 c5 c6 imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9 ...
    scenario cascade_scale_time age_mix start_year treat ...
    APRI


Y=reshape(y,3,10,9,20);

S=Y(:,:,:,1);
S1=Y(:,:,:,2);
S2=Y(:,:,:,3);
S3=Y(:,:,:,4);
S4=Y(:,:,:,5);
A=Y(:,:,:,6);
T=Y(:,:,:,7);
T1=Y(:,:,:,8);
T2=Y(:,:,:,9);
T3=Y(:,:,:,10);
T4=Y(:,:,:,11);
F0=Y(:,:,:,12);
F1=Y(:,:,:,13);
F2=Y(:,:,:,14);
F3=Y(:,:,:,15);
F4=Y(:,:,:,16);
DC=Y(:,:,:,17);
HCC=Y(:,:,:,18);
LT=Y(:,:,:,19);
LT2=Y(:,:,:,20);

pop=sum(sum(sum(Y(:,:,:,1:20),2),3),4);

R_inc=1;

%Version that gets ~9% in F34 at 2015
%if TT<50 R_inc=1+ (1.5/50)*(TT); elseif TT>=50 && TT<55 R_inc=2.5-(1.5/5)*(TT-50); elseif TT>=55 R_inc=1; end
%if TT<25 R_inc=1; elseif TT>=25 && TT<30 R_inc=1-(0.27/5)*(TT-25); elseif TT>=30 R_inc=0.72; end
%if TT<50 R_inc=1+ (1.5/50)*(TT); elseif TT>=50 R_inc=2.5; end
%R_inc = R_inc / 2.5;
%% Sensitivity analysis of R_inc
%if TT<40 R_inc=1; elseif TT>=40 && TT<50 R_inc=1+(1.5/10)*(TT-40); elseif TT>=50 && TT<55 R_inc=2.5-(1.5/5)*(TT-50); elseif TT>=55 R_inc=1; end %Peak only in 90s
%if TT<50 R_inc=1+ (3/50)*(TT); elseif TT>=50 && TT<55 R_inc=4-(3/5)*(TT-50); elseif TT>=55 R_inc=1; end %Double the height
%R_inc=1; %Flat
%if TT<50 R_inc=1.8; elseif TT>=50 && TT<55 R_inc=1.8-(0.5/5)*(TT-50); elseif TT>=55 R_inc=1.3; end %Constant then decrease
%if TT<40 R_inc=0.1+(0.9/40)*TT; elseif TT>=40 && TT<50 R_inc=1+(1.5/10)*(TT-40); elseif TT>=50 && TT<55 R_inc=2-(1.5/5)*(TT-50); elseif TT>=55 R_inc=1; end

%Plausible based on 250% increase in HCV notifications in the 1990s. 7% at
%F34 condition
%if TT<40 R_inc=1; elseif TT>=40 && TT<50 R_inc=1+(1.5/10)*(TT-40); elseif TT>=50 && TT<55 R_inc=2-(1.5/5)*(TT-50); elseif TT>=55 R_inc=1; end

%R_inc=rel_inc(find([5:5:60,1000]>=TT,1));

%%

    function ad = rate(r_XY_PWID,r_XY) % sub-function that distributes factors for PWID and former PWID
        ad = reshape(repmat([r_XY_PWID,r_XY,r_XY],9,10)',3,10,9);
    end

if TT > start_year
    % Force of infection
    lambda=R_inc*infect*sum(sum(A(1,:,:)+F0(1,:,:)+F1(1,:,:)+F2(1,:,:)+F3(1,:,:)+F4(1,:,:)+DC(1,:,:)+HCC(1,:,:)+LT(1,:,:)+LT2(1,:,:)))/sum(sum(sum(Y(1,:,:,1:20))));
    if age_mix < 0
        lambda2 = lambda;
    else
        lambda_age = zeros(1,9);
        for a=2:8
            if sum(sum(sum(Y(1,:,a,1:20)))) > 0
                lambda_age(a) = R_inc*infect*(...
                    age_mix * sum(sum(A(1,:,a)+F0(1,:,a)+F1(1,:,a)+F2(1,:,a)+F3(1,:,a)+F4(1,:,a)+DC(1,:,a)+HCC(1,:,a)+LT(1,:,a)+LT2(1,:,a)))/sum(sum(sum(Y(1,:,a,1:20)))) + ...
                    (1-age_mix)*sum(sum(A(1,:,[1:a-1,a+1:end])+F0(1,:,[1:a-1,a+1:end])+F1(1,:,[1:a-1,a+1:end])+F2(1,:,[1:a-1,a+1:end])+F3(1,:,[1:a-1,a+1:end])+F4(1,:,[1:a-1,a+1:end])+DC(1,:,[1:a-1,a+1:end])+HCC(1,:,[1:a-1,a+1:end])+LT(1,:,[1:a-1,a+1:end])+LT2(1,:,[1:a-1,a+1:end])))/sum(sum(sum(Y(1,:,[1:a-1,a+1:end],1:20)))) ...
                    );
            end
        end
        if sum(sum(sum(Y(1,:,1,1:20)))) > 0 && sum(sum(sum(Y(1,:,2:end,1:20)))) > 0
            lambda_age(1) = R_inc*infect*(...
                age_mix*sum(sum(A(1,:,1)+F0(1,:,1)+F1(1,:,1)+F2(1,:,1)+F3(1,:,1)+F4(1,:,1)+DC(1,:,1)+HCC(1,:,1)+LT(1,:,1)+LT2(1,:,1)))/sum(sum(sum(Y(1,:,1,1:20)))) + ...
                (1-age_mix)*sum(sum(A(1,:,2:end)+F0(1,:,2:end)+F1(1,:,2:end)+F2(1,:,2:end)+F3(1,:,2:end)+F4(1,:,2:end)+DC(1,:,2:end)+HCC(1,:,2:end)+LT(1,:,2:end)+LT2(1,:,2:end)))/sum(sum(sum(Y(1,:,2:end,1:20)))) ...
                );
        end
        if sum(sum(sum(Y(1,:,9,1:20)))) > 0
            lambda_age(9) = R_inc*infect*(...
                age_mix*sum(sum(A(1,:,9)+F0(1,:,9)+F1(1,:,9)+F2(1,:,9)+F3(1,:,9)+F4(1,:,9)+DC(1,:,9)+HCC(1,:,9)+LT(1,:,9)+LT2(1,:,9)))/sum(sum(sum(Y(1,:,9,1:20)))) + ...
                (1-age_mix)*sum(sum(A(1,:,1:end-1)+F0(1,:,1:end-1)+F1(1,:,1:end-1)+F2(1,:,1:end-1)+F3(1,:,1:end-1)+F4(1,:,1:end-1)+DC(1,:,1:end-1)+HCC(1,:,1:end-1)+LT(1,:,2:end)+LT2(1,:,1:end-1)))/sum(sum(sum(Y(1,:,1:end-1,1:20)))) ...
                );
        end
        lambda2 = cat(3,repmat([lambda_age(1);0;0],1,10),repmat([lambda_age(2);0;0],1,10),repmat([lambda_age(3);0;0],1,10),repmat([lambda_age(4);0;0],1,10),...
            repmat([lambda_age(5);0;0],1,10),repmat([lambda_age(6);0;0],1,10),repmat([lambda_age(7);0;0],1,10),repmat([lambda_age(8);0;0],1,10),repmat([lambda_age(9);0;0],1,10));
    end
    %         if TT>=75
    %             r_S4death = 0.002;
    %         end
    %TREATMENT ALLOCATION FUNCTION, f
    
    f=100*ones(3,10,9,20);
    if alpha < 0.8 && TT > 30
        %trea = sum(sum(A(1,:,:)+F0(1,:,:)+F1(1,:,:)+F2(1,:,:)+F3(1,:,:)+F4(1,:,:)+DC(1,:,:)+HCC(1,:,:)+LT(1,:,:)+LT2(1,:,:)));
        c5_PWID = -log(1-0.02*sum(sum(sum(Y(1:3,1,:,12:20))))/sum(sum(sum(Y(1:3,1,[6,8],1:20))))); %* sum(sum(A(1,:,:)+F0(1,:,:)+F1(1,:,:)+F2(1,:,:)+F3(1,:,:)+F4(1,:,:)+DC(1,:,:)+HCC(1,:,:)+LT(1,:,:)+LT2(1,:,:)));
        c5_formerPWID = -log(1-0.02*sum(sum(sum(Y(1:3,1,:,12:20))))/sum(sum(sum(Y(1:3,1,[6,8],1:20))))); %* sum(sum(A(2,:,:)+F0(2,:,:)+F1(2,:,:)+F2(2,:,:)+F3(2,:,:)+F4(2,:,:)+DC(2,:,:)+HCC(2,:,:)+LT(2,:,:)+LT2(2,:,:)));
        c5 = -log(1-0.02*sum(sum(sum(Y(1:3,1,:,12:20))))/sum(sum(sum(Y(1:3,1,[6,8],1:20))))); %* sum(sum(A(3,:,:)+F0(3,:,:)+F1(3,:,:)+F2(3,:,:)+F3(3,:,:)+F4(3,:,:)+DC(3,:,:)+HCC(3,:,:)+LT(3,:,:)+LT2(3,:,:)));
        c6_PWID = 0;%-log(1-0.02); %* sum(sum(A(1,:,:)+F0(1,:,:)+F1(1,:,:)+F2(1,:,:)+F3(1,:,:)+F4(1,:,:)+DC(1,:,:)+HCC(1,:,:)+LT(1,:,:)+LT2(1,:,:)));
        c6_formerPWID = 0;%-log(1-0.02); %* sum(sum(A(2,:,:)+F0(2,:,:)+F1(2,:,:)+F2(2,:,:)+F3(2,:,:)+F4(2,:,:)+DC(2,:,:)+HCC(2,:,:)+LT(2,:,:)+LT2(2,:,:)));
        c6 = 0;%-log(1-0.02); %* sum(sum(A(3,:,:)+F0(3,:,:)+F1(3,:,:)+F2(3,:,:)+F3(3,:,:)+F4(3,:,:)+DC(3,:,:)+HCC(3,:,:)+LT(3,:,:)+LT2(3,:,:)));
    end
    if TT>50 && TT<=66
        trea=treat(3);
        trea_extra=0;
        [f]=treatment_alloc(Y,trea,trea_extra);
    end
    if TT>66
        trea=treat(1);
        if TT<=67 trea_extra=treat(2); else trea_extra=0; end
        [f]=treatment_alloc(Y,trea,trea_extra);
    end
    
    new0 = lambda2.*S;
    new1 = lambda2.*S1;
    new2 = lambda2.*S2;
    new3 = lambda2.*S3;
    new4 = lambda2.*S4;
    
    %Infection and progression
    ydot1 = cat(4,delta*r_AF0*A-new0.*rate(1,0),... %S - to get infected must be PWID
        zeros(3,10,9),... %S1
        zeros(3,10,9),... %S2
        zeros(3,10,9),... %S3
        -r_S4death*S4,... %S4
        -r_AF0*A,... %Acute
        zeros(3,10,9),...%T
        zeros(3,10,9),...%T1
        zeros(3,10,9),...%T2
        zeros(3,10,9),...%T3
        -r_S4death*T4,...%T4
        (1-delta)*r_AF0*A-rate(r_F0F1_PWID,r_F0F1).*F0,...%F0
        rate(r_F0F1_PWID,r_F0F1).*F0-rate(r_F1F2_PWID,r_F1F2).*F1,...%F1
        rate(r_F1F2_PWID,r_F1F2).*F1-rate(r_F2F3_PWID,r_F2F3).*F2,...%F2
        rate(r_F2F3_PWID,r_F2F3).*F2-rate(r_F3F4_PWID,r_F3F4).*F3,...%F3
        rate(r_F3F4_PWID,r_F3F4).*F3-r_F4DC*F4-r_F4HCC*F4,...%F4
        r_F4DC*F4-r_DCHCC*DC-r_DCLT*DC-r_DCdeath*DC,...%DC
        r_DCHCC*DC+r_F4HCC*F4-r_HCCLT*HCC-r_HCCdeath*HCC,...%HCC
        r_HCCLT*HCC+r_DCLT*DC-r_LTdeath1*LT-r_LT1LT2*LT,... %LT
        r_LT1LT2*LT-r_LTdeath2*LT2); %LT2
    
    
    
    %Mortality
    MU = reshape(repmat([mu_PWID;mu_former;mu_former],10,1),3,10,9).*rate(R_inc,1); % assumes current and former die at same rate
    l_death=sum(sum(sum(r_HCCdeath*HCC(:,:,:)+r_LTdeath1*LT(:,:,:)+r_LTdeath2*LT2(:,:,:)+r_DCdeath*DC(:,:,:)+r_S4death*S4(:,:,:)+r_S4death*T4(:,:,:))));
    ydot3=cat(4,-S.*MU,...
        -S1.*MU,-S2.*MU,-S3.*MU,-S4.*MU,-A.*MU,-T.*MU,-T1.*MU,-T2.*MU,-T3.*MU,-T4.*MU,-F0.*MU,-F1.*MU,-F2.*MU,-F3.*MU,-F4.*MU,-DC.*MU,-HCC.*MU,-LT.*MU,-LT2.*MU);
    
    
    %Cessation and relapse into IDU
    r_relapse_rel=r_relapse;%*R_inc;
    exit_IDU_rel=exit_IDU;%/R_inc;
    
    ydot4=zeros(3,10,9,20);
    ydot4(1,:,:,1:20)=-Y(1,:,:,1:20)*exit_IDU_rel + Y(2,:,:,1:20)*r_relapse;
    ydot4(2,:,:,1:20)=Y(1,:,:,1:20)*exit_IDU_rel-Y(2,:,:,1:20)*r_relapse- 1/(10*12)*Y(2,:,:,1:20); % former for 10 years
    ydot4(3,:,:,1:20)=1/(10*12)*Y(2,:,:,1:20);
    
    
    %Let people age
    ydot5=zeros(3,10,9,20);
    ydot5(:,:,1,1:20)=-1/5*Y(:,:,1,1:20);
    ydot5(:,:,2,1:20)=1/5*Y(:,:,1,1:20)-1/5*Y(:,:,2,1:20);
    ydot5(:,:,3,1:20)=1/5*Y(:,:,2,1:20)-1/5*Y(:,:,3,1:20);
    ydot5(:,:,4,1:20)=1/5*Y(:,:,3,1:20)-1/10*Y(:,:,4,1:20);
    ydot5(:,:,5,1:20)=1/10*Y(:,:,4,1:20)-1/10*Y(:,:,5,1:20);
    ydot5(:,:,6,1:20)=1/10*Y(:,:,5,1:20)-1/10*Y(:,:,6,1:20);
    ydot5(:,:,7,1:20)=1/10*Y(:,:,6,1:20)-1/10*Y(:,:,7,1:20);
    ydot5(:,:,8,1:20)=1/10*Y(:,:,7,1:20)-1/10*Y(:,:,8,1:20);
    ydot5(:,:,9,1:20)=1/10*Y(:,:,8,1:20);
    
    %New arrivals
    IDUs=PWID0;
    if TT<50
        IDUs = (TT-start_year)*(400-PWID0)/(50-start_year) + PWID0;
    elseif TT>=50 && TT<55
        IDUs = 400 - (TT-50)*(400 -total_PWID)/5;
    elseif TT>=55
        IDUs = total_PWID;
    end
    mu_arrivals = max(IDUs-sum(sum(sum(Y(1,:,:,1:20))))...
        +sum(sum(sum(-ydot3(1,:,:,:))))+sum(sum(sum(-ydot4(1,:,:,1:20)))),...
        0);
    %         mu_arrivals=sum(sum(sum(sum(-ydot3(:,:,:,:)))));
    ydot6=zeros(3,10,9,20);
    ydot6(1,1,1,1)=l_death+mu_arrivals;
    %             import_infections=0;
    %             if TT<25 import_infections=imp1;
    %             elseif TT>=25 && TT<30 import_infections=imp2;
    %             elseif TT>=30 && TT<35 import_infections=imp3;
    %             elseif TT>=35 && TT<40 import_infections=imp4;
    %             elseif TT>=40 && TT<45 import_infections=imp5;
    %             elseif TT>=45 && TT<50 import_infections=imp6;
    %             elseif TT>=50 && TT<55 import_infections=imp7;
    %             elseif TT>=55 && TT<60 import_infections=imp8;
    %             elseif TT>=60 && TT<65 import_infections=imp9;
    %             elseif TT>=65 import_infections=imp9; end
    %             ydot6(2,1,1,12)=0.5*import_infections;
    %             ydot6(2,1,1,1)=(0.5*import_infections)/max(1,(sum(sum(sum(Y(2,:,:,12:20))))/max(1,sum(sum(sum(Y(2,:,:,1:20))))))); % keep prevalence constant among former
    %             ydot6(3,1,1,12)=0.5*import_infections;
    %
    %             ydot1(2,1,1,27) = 0.5*import_infections;
    %             ydot1(3,1,1,27) = 0.5*import_infections;
    
    %% Progression in cascade
    if TT > 66 && TT <= 67
        mult = min(1, (TT-66));
        apply_cascade_rates('initial', mult);
    end
    if TT > 67 && TT <= 67.5 + cascade_scale_time
        if cascade_scale_time == 0
            mult = 1;
        else
            mult = min(1, 1 - ((cascade_scale_time - (TT-67)) / cascade_scale_time));
        end
        apply_cascade_rates(scenario, mult);
    end
    
    %Current PWID
    ydot7=zeros(3,10,9,20);
    
    ydot7(1,1:4,:,17:20) = -12*Y(1,1:4,:,17:20); % assumes DC, HCC and LT HCV antibodies are detected and made ready for treatment immediately
    ydot7(1,6,:,17:20) = reshape(sum(12*Y(1,1:4,:,17:20),2),1,1,9,4);
    
    ydot7(1,1,:,12:16) = -c1_chronicPWID*Y(1,1,:,12:16); % HCV antibody screening rate for F0-F4 disease
    ydot7(1,2,:,12:16) = c1_chronicPWID*Y(1,1,:,12:16) - c2_PWID*Y(1,2,:,12:16); % HCV RNA+
    ydot7(1,3,:,12:16) = c2_PWID*Y(1,2,:,12:16) - c3_PWID*Y(1,3,:,12:16); % PCR test for genotype
    ydot7(1,4,:,12:16) = c3_PWID*Y(1,3,:,12:16) - c4_PWID*Y(1,4,:,12:16); % Liver disease test (Fibroscan)
    ydot7(1,6,:,12:16) = c4_PWID*Y(1,4,:,12:16); % ready for treatment

    
    %Former PWID
    ydot7(2,1:4,:,17:20) = -12*Y(2,1:4,:,17:20); % assumes DC, HCC and LT HCV antibodies are detected and made ready for treatment immediately
    ydot7(2,6,:,17:20) = reshape(sum(12*Y(2,1:4,:,17:20),2),1,1,9,4);
    
    ydot7(2,1,:,12:16) = -c1_chronicformerPWID*Y(2,1,:,12:16); % HCV antibody screening rate for F0-F4 disease
    ydot7(2,2,:,12:16) = c1_chronicformerPWID*Y(2,1,:,12:16) - c2_formerPWID*Y(2,2,:,12:16); % HCV RNA+
    ydot7(2,3,:,12:16) = c2_formerPWID*Y(2,2,:,12:16) - c3_formerPWID*Y(2,3,:,12:16); % PCR test for genotype
    ydot7(2,4,:,12:16) = c3_formerPWID*Y(2,3,:,12:16) - c4_formerPWID*Y(2,4,:,12:16); % Liver disease test (Fibroscan)
    ydot7(2,6,:,12:16) = c4_formerPWID*Y(2,4,:,12:16); % ready for treatment
        
    %Never injected
    ydot7(3,1:4,:,17:20) = -12*Y(3,1:4,:,17:20); % assumes DC, HCC and LT HCV antibodies are detected and made ready for treatment immediately
    ydot7(3,6,:,17:20) = reshape(sum(12*Y(3,1:4,:,17:20),2),1,1,9,4);
    
    ydot7(3,1,:,12:16) = -c1_chronic*Y(3,1,:,12:16); % HCV antibody screening rate for F0-F4 disease
    ydot7(3,2,:,12:16) = c1_chronic*Y(3,1,:,12:16) - c2*Y(3,2,:,12:16); % HCV RNA+
    ydot7(3,3,:,12:16) = c2*Y(3,2,:,12:16) - c3*Y(3,3,:,12:16); % test for genotype
    ydot7(3,4,:,12:16) = c3*Y(3,3,:,12:16) - c4*Y(3,4,:,12:16); % Liver disease test (Fibroscan)
    ydot7(3,6,:,12:16) = c4*Y(3,4,:,12:16); % ready for treatment

    
    
    if TT > 67 && APRI == 1 % Fibroscan not required for F2 or less if scenario uses APRI
        if cascade_scale_time==0 mult=1; else mult = min(1, 1 - ((cascade_scale_time - (TT-67)) / cascade_scale_time)); end
        ydot7(1,4,:,12:14) = c3_PWID*Y(1,3,:,12:14) - ((1-mult)*c4_PWID+mult*52)*Y(1,4,:,12:14); % Fibroscan skipped for F0-F2
        ydot7(1,6,:,12:14) = ((1-mult)*c4_PWID+mult*52)*Y(1,4,:,12:14); % ready for treatment
        ydot7(2,4,:,12:14) = c3_formerPWID*Y(2,3,:,12:14) - ((1-mult)*c4_formerPWID+mult*52)*Y(2,4,:,12:14); % Fibroscan skipped for F0-F2
        ydot7(2,6,:,12:14) = ((1-mult)*c4_formerPWID+mult*52)*Y(2,4,:,12:14); % ready for treatment
        ydot7(3,4,:,12:14) = c3*Y(3,3,:,12:14) - ((1-mult)*c4+mult*52)*Y(3,4,:,12:14); % Fibroscan skipped for F0-F2
        ydot7(3,6,:,12:14) = ((1-mult)*c4+mult*52)*Y(3,4,:,12:14); % ready for treatment

    end
    
    %Treatment
    ydot8=zeros(3,10,9,20);
    ydot8(:,6,:,12) = - min(f(:,6,:,12), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*F0(:,6,:)); % treatment from F0
    ydot8(:,6,:,13) = - min(f(:,6,:,13), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*F1(:,6,:)); % treatment from F1
    ydot8(:,6,:,14) = - min(f(:,6,:,14), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*F2(:,6,:)); % treatment from F2
    ydot8(:,6,:,15) = - min(f(:,6,:,15), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*F3(:,6,:)); % treatment from F3
    ydot8(:,6,:,16) = - min(f(:,6,:,16), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*F4(:,6,:)); % treatment from F4
    ydot8(:,6,:,17) = - min(f(:,6,:,17), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*DC(:,6,:)); % treatment from DC
    ydot8(:,6,:,18) = - min(f(:,6,:,18), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*HCC(:,6,:)); % treatment from HCC
    ydot8(:,6,:,19) = - min(f(:,6,:,19), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*LT(:,6,:)); % treatment from LT
    ydot8(:,6,:,20) = - min(f(:,6,:,20), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*LT2(:,6,:)); % treatment from LT2
    ydot8(:,7,:,7) = reshape(repmat([alpha*p_complete,alpha,alpha],9,1)',3,1,9).*(-ydot8(:,6,:,12)); % in treatment from F0 and going to succeed
    ydot8(:,7,:,8) = reshape(repmat([alpha*p_complete,alpha,alpha],9,1)',3,1,9).*(-ydot8(:,6,:,13)); % T1
    ydot8(:,7,:,9) = reshape(repmat([alpha*p_complete,alpha,alpha],9,1)',3,1,9).*(-ydot8(:,6,:,14)); % T2
    ydot8(:,7,:,10) = reshape(repmat([alpha*p_complete,alpha,alpha],9,1)',3,1,9).*(-ydot8(:,6,:,15)); % T3
    ydot8(:,7,:,11) = reshape(repmat([alpha*p_complete,alpha,alpha],9,1)',3,1,9).*sum((-ydot8(:,6,:,16:20)),4); % T4
    ydot8(:,8,:,12) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot8(:,6,:,12)); % failing treatment from F0
    ydot8(:,8,:,13) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot8(:,6,:,13)); % failing treatment from F1
    ydot8(:,8,:,14) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot8(:,6,:,14)); % failing treatment from F2
    ydot8(:,8,:,15) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot8(:,6,:,15)); % failing treatment from F3
    ydot8(:,8,:,16) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot8(:,6,:,16)); % failing treatment from F4
    ydot8(:,8,:,17) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot8(:,6,:,17));% failing treatment from DC
    ydot8(:,8,:,18) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot8(:,6,:,18));% failing treatment from HCC
    ydot8(:,8,:,19) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot8(:,6,:,19)); % failing treatment from LT
    ydot8(:,8,:,20) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot8(:,6,:,20)); % failing treatment from LT2
    ydot8(:,7,:,7) = ydot8(:,7,:,7) - omega*T(:,7,:); % Recovering and moving to the 1-st level
    ydot8(:,7,:,8) = ydot8(:,7,:,8) - omega*T1(:,7,:); % Recovering and moving to the 1-st level
    ydot8(:,7,:,9) = ydot8(:,7,:,9) - omega*T2(:,7,:); % Recovering and moving to the 1-st level
    ydot8(:,7,:,10) = ydot8(:,7,:,10) - omega*T3(:,7,:); % Recovering and moving to the 1-st level
    ydot8(:,7,:,11) = ydot8(:,7,:,11) - omega*T4(:,7,:); % Recovering and moving to the 1-st level
    ydot8(:,1,:,1) = omega*T(:,7,:); % Recovering and moving to the 1-st level
    ydot8(:,1,:,2) = omega*T1(:,7,:); % Recovering and moving to the 1-st level
    ydot8(:,1,:,3) = omega*T2(:,7,:); % Recovering and moving to the 1-st level
    ydot8(:,1,:,4) = omega*T3(:,7,:); % Recovering and moving to the 1-st level
    ydot8(:,1,:,5) = omega*T4(:,7,:); % Recovering and moving to the 1-st level

    
    
    %Retreatment
    ydot9=zeros(3,10,9,20);
    ydot9(:,8,:,12) = - min(f(:,8,:,12), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*F0(:,8,:)); % treatment from F0
    ydot9(:,8,:,13) = - min(f(:,8,:,13), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*F1(:,8,:)); % treatment from F1
    ydot9(:,8,:,14) = - min(f(:,8,:,14), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*F2(:,8,:)); % treatment from F2
    ydot9(:,8,:,15) = - min(f(:,8,:,15), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*F3(:,8,:)); % treatment from F3
    ydot9(:,8,:,16) = - min(f(:,8,:,16), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*F4(:,8,:)); % treatment from F4
    ydot9(:,8,:,17) = - min(f(:,8,:,17), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*DC(:,8,:)); % treatment from DC
    ydot9(:,8,:,18) = - min(f(:,8,:,18), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*HCC(:,8,:)); % treatment from HCC
    ydot9(:,8,:,19) = - min(f(:,8,:,19), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*LT(:,8,:)); % treatment from LT
    ydot9(:,8,:,20) = - min(f(:,8,:,20), reshape(repmat([c5_PWID,c5_formerPWID,c5],9,1)',3,1,9).*LT2(:,8,:)); % treatment from LT2
    ydot9(:,9,:,7) = reshape(repmat([alpha*p_complete,alpha,alpha],9,1)',3,1,9).*(-ydot9(:,8,:,12)); % in treatment from F0 and going to succeed
    ydot9(:,9,:,8) = reshape(repmat([alpha*p_complete,alpha,alpha],9,1)',3,1,9).*(-ydot9(:,8,:,13)); % T1
    ydot9(:,9,:,9) = reshape(repmat([alpha*p_complete,alpha,alpha],9,1)',3,1,9).*(-ydot9(:,8,:,14)); % T2
    ydot9(:,9,:,10) = reshape(repmat([alpha*p_complete,alpha,alpha],9,1)',3,1,9).*(-ydot9(:,8,:,15)); % T3
    ydot9(:,9,:,11) = reshape(repmat([alpha*p_complete,alpha,alpha],9,1)',3,1,9).*sum((-ydot9(:,8,:,16:20)),4); % T4
    ydot9(:,10,:,12) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot9(:,8,:,12)); % failing treatment from F0
    ydot9(:,10,:,13) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot9(:,8,:,13)); % failing treatment from F1
    ydot9(:,10,:,14) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot9(:,8,:,14)); % failing treatment from F2
    ydot9(:,10,:,15) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot9(:,8,:,15)); % failing treatment from F3
    ydot9(:,10,:,16) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot9(:,8,:,16)); % failing treatment from F4
    ydot9(:,10,:,17) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot9(:,8,:,17));% failing treatment from DC
    ydot9(:,10,:,18) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot9(:,8,:,18));% failing treatment from HCC
    ydot9(:,10,:,19) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot9(:,8,:,19)); % failing treatment from LT
    ydot9(:,10,:,20) = reshape(repmat([1-alpha*p_complete,1-alpha,1-alpha],9,1)',3,1,9).*(-ydot9(:,8,:,20)); % failing treatment from LT2
    ydot9(:,9,:,7) = ydot9(:,9,:,7) - omega*T(:,9,:); % Recovering and moving to the 1-st level
    ydot9(:,9,:,8) = ydot9(:,9,:,8) - omega*T1(:,9,:); % Recovering and moving to the 1-st level
    ydot9(:,9,:,9) = ydot9(:,9,:,9) - omega*T2(:,9,:); % Recovering and moving to the 1-st level
    ydot9(:,9,:,10) = ydot9(:,9,:,10) - omega*T3(:,9,:); % Recovering and moving to the 1-st level
    ydot9(:,9,:,11) = ydot9(:,9,:,11) - omega*T4(:,9,:); % Recovering and moving to the 1-st level
    ydot9(:,1,:,1) = omega*T(:,9,:); % Recovering and moving to the 1-st level
    ydot9(:,1,:,2) = omega*T1(:,9,:); % Recovering and moving to the 1-st level
    ydot9(:,1,:,3) = omega*T2(:,9,:); % Recovering and moving to the 1-st level
    ydot9(:,1,:,4) = omega*T3(:,9,:); % Recovering and moving to the 1-st level
    ydot9(:,1,:,5) = omega*T4(:,9,:); % Recovering and moving to the 1-st level

    
    
    ydot=reshape(ydot1+ydot3+ydot4+ydot5+ydot6+ydot7+ydot8+ydot9,3*10*9*20,1);
    
else
    ydot = zeros(3*10*9*20,1);
end

end