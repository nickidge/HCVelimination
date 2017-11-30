
function[TT,y]=DE_track_age(Tin,y0,t0,treat)

global mu_PWID mu_former exit_IDU r_relapse delta alpha p_complete omega infect total_PWID PWID0 dem dt...
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2...
    imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9 imported...
    scenario cascade_scale_time age_mix start_year r_inc_up followup ...
    APRI num_pops num_cascade num_age num_intervention num_engagement num_region infect_factor progression progression_base...
    ost_enrollment ost_duration nsp_enrollment nsp_duration RNAtesting harm_reduction_coverage ost_coverage nsp_coverage


APRI = 1;

y0=reshape(y0,num_pops*num_cascade*num_age*num_intervention*num_engagement*num_region*(27+6),1); % adding cascade states to capture transfers

tvec = [t0(end):dt:t0(end)+Tin]'; % time steps

y(1,:) = [y0];
%ODE subfunction
    for t = 2:length(tvec)
        TT = tvec(t);
        if TT > start_year
            
            Y=reshape(y(t-1,:),num_pops, num_cascade, num_age, num_intervention, num_engagement,num_region,27+6);
            
            S=Y(:,:,:,:,:,:,1);
            S1=Y(:,:,:,:,:,:,2);
            S2=Y(:,:,:,:,:,:,3);
            S3=Y(:,:,:,:,:,:,4);
            S4=Y(:,:,:,:,:,:,5);
            A=Y(:,:,:,:,:,:,6);
            T=Y(:,:,:,:,:,:,7);
            T1=Y(:,:,:,:,:,:,8);
            T2=Y(:,:,:,:,:,:,9);
            T3=Y(:,:,:,:,:,:,10);
            T4=Y(:,:,:,:,:,:,11);
            F0=Y(:,:,:,:,:,:,12);
            F1=Y(:,:,:,:,:,:,13);
            F2=Y(:,:,:,:,:,:,14);
            F3=Y(:,:,:,:,:,:,15);
            F4=Y(:,:,:,:,:,:,16);
            DC=Y(:,:,:,:,:,:,17);
            HCC=Y(:,:,:,:,:,:,18);
            LT=Y(:,:,:,:,:,:,19);
            LT2=Y(:,:,:,:,:,:,20);
            F4_transfer=Y(:,:,:,:,:,:,21);
            Ldeath=Y(:,:,:,:,:,:,22);
            T_total=Y(:,:,:,:,:,:,23);
            HCC_transfer=Y(:,:,:,:,:,:,24);
            T_F4on_total=Y(:,:,:,:,:,:,25);
            Liver_transplants=Y(:,:,:,:,:,:,26);
            Inc=Y(:,:,:,:,:,:,27);
            
            Cas1=Y(:,:,:,:,:,:,28);
            Cas2=Y(:,:,:,:,:,:,29);
            Cas3=Y(:,:,:,:,:,:,30);
            Cas4=Y(:,:,:,:,:,:,31);
            Cas5=Y(:,:,:,:,:,:,32);
            Cas6=Y(:,:,:,:,:,:,33);
            
            pop=reshape(sum(sum(sum(sum(sum(Y(:,:,:,:,:,:,1:20),2),3),4),5),7),num_pops, num_region); %vector of population sizes
            
            %% Relative incidence function
            R_inc=1;
            
            %Version that gets ~9% in F34 at 2015
            if TT<40 R_inc=1+ (r_inc_up/40)*(TT); elseif TT>=40 && TT<45 R_inc=1+r_inc_up-(r_inc_up/5)*(TT-40); elseif TT>=45 R_inc=1; end
            
            
            
            %% Force of infection
            
            lambda_age = zeros(num_pops,num_age,num_region);
            for j = 1:num_region
                for i =1:num_pops
                    %                        if i==1 % identifying which population groups mix
                    %                             p = [1,5]; % all PWID population groups
                    %                         elseif i==4 || i==6
                    %                             p = [4,5,6]; % all MSM population groups
                    %                         elseif i ==5 % MSM and PWID
                    %                             p = [1,4,5,6];
                    %                         else
                    p = i;
                    %                        end
                    if i == 4
                        for a = 1:num_age
                            lambda_age(i,a,j) = max(0,infect_MSM*(...
                            sum(sum(sum(sum(A(p,:,:,:,:,j)+F0(p,:,:,:,:,j)+F1(p,:,:,:,:,j)+F2(p,:,:,:,:,j)+F3(p,:,:,:,:,j)+F4(p,:,:,:,:,j)+DC(p,:,:,:,:,j)+HCC(p,:,:,:,:,j)+LT(p,:,:,:,:,j)+LT2(p,:,:,:,:,j)))))/sum(sum(sum(sum(sum(Y(p,:,:,:,:,j,1:20)))))) ...
                            ));
                        end
                    end
                    if sum(sum(sum(sum(Y(i,:,1,:,:,j,1:20))))) > 0 && sum(sum(sum(sum(sum(Y(i,:,2:end,:,:,j,1:20)))))) > 0
                        lambda_age(i,1,j) = max(0,R_inc*infect*(...
                            age_mix*sum(sum(sum(sum(A(p,:,1,:,:,j)+F0(p,:,1,:,:,j)+F1(p,:,1,:,:,j)+F2(p,:,1,:,:,j)+F3(p,:,1,:,:,j)+F4(p,:,1,:,:,j)+DC(p,:,1,:,:,j)+HCC(p,:,1,:,:,j)+LT(p,:,1,:,:,j)+LT2(p,:,1,:,:,j)))))/sum(sum(sum(sum(sum(Y(p,:,1,:,:,j,1:20)))))) + ...
                            (1-age_mix)*sum(sum(sum(sum(sum(A(p,:,2:end,:,:,j)+F0(p,:,2:end,:,:,j)+F1(p,:,2:end,:,:,j)+F2(p,:,2:end,:,:,j)+F3(p,:,2:end,:,:,j)+F4(p,:,2:end,:,:,j)+DC(p,:,2:end,:,:,j)+HCC(p,:,2:end,:,:,j)+LT(p,:,2:end,:,:,j)+LT2(p,:,2:end,:,:,j))))))/sum(sum(sum(sum(sum(sum(Y(p,:,2:end,:,:,j,1:20))))))) ...
                            ));
                    end
                    for a=2:num_age-1
                        if sum(sum(sum(sum(sum(Y(i,:,a,:,:,:,1:20)))))) > 0
                            lambda_age(i,a,j) = max(0, R_inc*infect*(...
                                age_mix * sum(sum(sum(sum(A(p,:,a,:,:,j)+F0(p,:,a,:,:,j)+F1(p,:,a,:,:,j)+F2(p,:,a,:,:,j)+F3(p,:,a,:,:,j)+F4(p,:,a,:,:,j)+DC(p,:,a,:,:,j)+HCC(p,:,a,:,:,j)+LT(p,:,a,:,:,j)+LT2(p,:,a,:,:,j)))))/sum(sum(sum(sum(sum(Y(p,:,a,:,:,j,1:20)))))) + ...
                                (1-age_mix)*sum(sum(sum(sum(sum(A(p,:,[1:a-1,a+1:end],:,:,j)+F0(p,:,[1:a-1,a+1:end],:,:,j)+F1(p,:,[1:a-1,a+1:end],:,:,j)+F2(p,:,[1:a-1,a+1:end],:,:,j)+F3(p,:,[1:a-1,a+1:end],:,:,j)+F4(p,:,[1:a-1,a+1:end],:,:,j)+DC(p,:,[1:a-1,a+1:end],:,:,j)+HCC(p,:,[1:a-1,a+1:end],:,:,j)+LT(p,:,[1:a-1,a+1:end],:,:,j)+LT2(p,:,[1:a-1,a+1:end],:,:,j))))))/sum(sum(sum(sum(sum(sum(Y(p,:,[1:a-1,a+1:end],:,:,j,1:20))))))) ...
                                ));
                        end
                    end
                    if sum(sum(sum(Y(i,:,num_age,1:20)))) > 0
                        lambda_age(i,num_age,j) = max(0, R_inc*infect*(...
                            age_mix*sum(sum(sum(sum(A(p,:,num_age,:,:,j)+F0(p,:,num_age,:,:,j)+F1(p,:,num_age,:,:,j)+F2(p,:,num_age,:,:,j)+F3(p,:,num_age,:,:,j)+F4(p,:,num_age,:,:,j)+DC(p,:,num_age,:,:,j)+HCC(p,:,num_age,:,:,j)+LT(p,:,num_age,:,:,j)+LT2(p,:,num_age,:,:,j)))))/sum(sum(sum(sum(sum(Y(p,:,num_age,:,:,j,1:20)))))) + ...
                            (1-age_mix)*sum(sum(sum(sum(sum(A(p,:,1:end-1,:,:,j)+F0(p,:,1:end-1,:,:,j)+F1(p,:,1:end-1,:,:,j)+F2(p,:,1:end-1,:,:,j)+F3(p,:,1:end-1,:,:,j)+F4(p,:,1:end-1,:,:,j)+DC(p,:,1:end-1,:,:,j)+HCC(p,:,1:end-1,:,:,j)+LT(p,:,1:end-1,:,:,j)+LT2(p,:,1:end-1,:,:,j))))))/sum(sum(sum(sum(sum(sum(Y(p,:,1:end-1,:,:,j,1:20))))))) ...
                            ));
                    end
                end
            end
            % reshape to match Y
            lambda2 = permute(reshape(repmat(infect_factor,1,num_age*num_cascade*num_engagement),num_pops,num_intervention, num_region,num_age, num_cascade,num_engagement),[1,5,4,2,6,3]) ...
                .* permute(reshape(repmat(lambda_age,1,num_cascade*num_engagement*num_intervention),num_pops,num_age,num_region,num_cascade,num_intervention,num_engagement),[1,4,2,5,6,3]);
            %lambda2 = lambda_age;
            
            %% Treatment allocation function, f
            
            f=100*ones(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region,20);
            if TT>52 && TT<=66 % TT = 52 corresponde to 2002 (1950+52); TT = 66 is 2016.
                trea=treat(3);
                trea_extra=0;
                [f]=treatment_alloc(Y,trea,trea_extra,0);
            end
            if TT>66
                trea=treat(1);
                if TT<=67 trea_extra=treat(2); else trea_extra=0; end
                [f]=treatment_alloc(Y,trea,trea_extra,0);
            end
            
            %% Infection and disease progression progression
            if TT>=66 imported2 = rate(imported,0); else imported2 = rate(0,0); end
            new0 = (lambda2).*S;
            new1 = (lambda2).*S1;
            new2 = (lambda2).*S2;
            new3 = (lambda2).*S3;
            new4 = (lambda2).*S4;
            
            ydot1 = 0*Y;
            ydot1 = cat(7,delta*r_AF0*A-new0,... %S
                -new1,... %S1
                -new2,... %S2
                -new3,... %S3
                -new4-r_S4death*S4,... %S4
                new0-r_AF0*A,... %Acute
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region),...%T
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region),...%T1
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region),...%T2
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region),...%T3
                -r_S4death*T4,...%T4
                (1-delta)*r_AF0*A-rate(r_F0F1_PWID,r_F0F1).*F0,...%F0
                new1+rate(r_F0F1_PWID,r_F0F1).*F0-rate(r_F1F2_PWID,r_F1F2).*F1,...%F1
                new2+rate(r_F1F2_PWID,r_F1F2).*F1-rate(r_F2F3_PWID,r_F2F3).*F2,...%F2
                new3+rate(r_F2F3_PWID,r_F2F3).*F2-rate(r_F3F4_PWID,r_F3F4).*F3,...%F3
                new4+rate(r_F3F4_PWID,r_F3F4).*F3-r_F4DC*F4-r_F4HCC*F4,...%F4
                r_F4DC*F4-r_DCHCC*DC-r_DCLT*DC-r_DCdeath*DC,...%DC
                r_DCHCC*DC+r_F4HCC*F4-r_HCCLT*HCC-r_HCCdeath*HCC,...%HCC
                r_HCCLT*HCC+r_DCLT*DC-r_LTdeath1*LT-r_LT1LT2*LT,... %LT
                r_LT1LT2*LT-r_LTdeath2*LT2,... %LT2
                new4+rate(r_F3F4_PWID,r_F3F4).*F3,... %Total transfers to F4 condition
                r_HCCdeath*HCC+r_LTdeath1*LT+r_LTdeath2*LT2+r_DCdeath*DC+r_S4death*S4+r_S4death*T4,...
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region),... %Total treatments from F0-F3
                r_DCHCC*DC+r_F4HCC*F4,... %Total transfers to HCC condition
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region),... %Total treatments from F4 onwards stage
                r_HCCLT*HCC+r_DCLT*DC,... %Total liver transplants
                new0+new1+new2+new3+new4,...%Total incidence
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region),... %Total movements c1-c2 (HCV antibody detection) in the cascade
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region),... %Total movements c2-c3 (HCV RNA detection) in the cascade
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region),... %Total movements c3-c4 (genotype) in the cascade
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region),... %Total movements c4-c5 (liver disease stage) in the cascade
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region),... %Total movements c5-c6
                zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region)); %Total movements c6-c7
            
            
            %% Mortality
            MU = permute(reshape(repmat([mu_PWID;mu_former;mu_former],1,num_cascade*num_intervention*num_engagement*num_region),...
                num_pops,num_age,num_cascade,num_intervention,num_engagement,num_region),[1,3,2,4,5,6]); % assumes current and former die at same rate
            l_death=sum(sum(sum(sum(sum(sum(r_HCCdeath*HCC+r_LTdeath1*LT+r_LTdeath2*LT2+r_DCdeath*DC+r_S4death*S4+r_S4death*T4))))));
            ydot3 = 0*Y;
            ydot3(:,:,:,:,:,:,1:20)=cat(7,-S.*MU,...
                -S1.*MU,-S2.*MU,-S3.*MU,-S4.*MU,-A.*MU,-T.*MU,-T1.*MU,-T2.*MU,-T3.*MU,-T4.*MU,...
                -F0.*MU,-F1.*MU,-F2.*MU,-F3.*MU,-F4.*MU,-DC.*MU,-HCC.*MU,-LT.*MU,-LT2.*MU);
            
            
            %% Cessation and relapse into IDU
            r_relapse_rel=r_relapse;%*R_inc;
            exit_IDU_rel=exit_IDU;%/R_inc;
            
            ydot4=0*Y;
            ydot4([1],:,:,:,:,:,1:20)=-Y([1],:,:,:,:,:,1:20)*exit_IDU_rel + Y([2],:,:,:,:,:,1:20)*r_relapse;
            ydot4([2],:,:,:,:,:,1:20)=Y([1],:,:,:,:,:,1:20)*exit_IDU_rel-Y([2],:,:,:,:,:,1:20)*r_relapse- 1/(5*12)*Y([2],:,:,:,:,:,1:20); % former for 10 years
            ydot4([3],:,:,:,:,:,1:20)=1/(5*12)*Y([2],:,:,:,:,:,1:20);
            
            
            %% Ageing
            ydot5=0*Y;
            ydot5(:,:,1,:,:,:,1:20)=-1/5*Y(:,:,1,:,:,:,1:20);
            ydot5(:,:,2,:,:,:,1:20)=1/5*Y(:,:,1,:,:,:,1:20)-1/5*Y(:,:,2,:,:,:,1:20);
            ydot5(:,:,3,:,:,:,1:20)=1/5*Y(:,:,2,:,:,:,1:20)-1/5*Y(:,:,3,:,:,:,1:20);
            ydot5(:,:,4,:,:,:,1:20)=1/5*Y(:,:,3,:,:,:,1:20)-1/10*Y(:,:,4,:,:,:,1:20);
            ydot5(:,:,5,:,:,:,1:20)=1/10*Y(:,:,4,:,:,:,1:20)-1/10*Y(:,:,5,:,:,:,1:20);
            ydot5(:,:,6,:,:,:,1:20)=1/10*Y(:,:,5,:,:,:,1:20)-1/10*Y(:,:,6,:,:,:,1:20);
            ydot5(:,:,7,:,:,:,1:20)=1/10*Y(:,:,6,:,:,:,1:20)-1/10*Y(:,:,7,:,:,:,1:20);
            ydot5(:,:,8,:,:,:,1:20)=1/10*Y(:,:,7,:,:,:,1:20)-1/10*Y(:,:,8,:,:,:,1:20);
            ydot5(:,:,9,:,:,:,1:20)=1/10*Y(:,:,8,:,:,:,1:20);
            
            %% New arrivals
            IDUs=PWID0;
            if TT<50 % TT = 50 is 2000
                IDUs = (TT-start_year)*(50000-PWID0)/(50-start_year) + PWID0;
            elseif TT>=50 && TT<55
                IDUs = 50000 - (TT-50)*(50000 -total_PWID)/5;
            elseif TT>=55
                IDUs = total_PWID;
            end
            mu_arrivals = max(IDUs-sum(sum(sum(sum(sum(sum(sum(Y([1],:,:,:,:,:,1:20))))))))...
                +sum(sum(sum(sum(sum(sum(sum(-ydot3([1],:,:,:,:,:,:))))))))+sum(sum(sum(sum(sum(sum(sum(-ydot4([1],:,:,:,:,:,1:20)))))))),...
                0);
            %         mu_arrivals=sum(sum(sum(sum(-ydot3(:,:,:,:)))));
            ydot6=0*Y;
            ydot6(1,1,1,1,2,1,1)=l_death+mu_arrivals;
            import_infections=0;
            if TT<25 import_infections=imp1;
            elseif TT>=25 && TT<30 import_infections=imp2;
            elseif TT>=30 && TT<35 import_infections=imp3;
            elseif TT>=35 && TT<40 import_infections=imp4;
            elseif TT>=40 && TT<45 import_infections=imp5;
            elseif TT>=45 && TT<50 import_infections=imp6;
            elseif TT>=50 && TT<55 import_infections=imp7;
            elseif TT>=55 && TT<60 import_infections=imp8;
            elseif TT>=60 && TT<65 import_infections=imp9;
            elseif TT>=65 import_infections=0; end
            ydot6(2,1,1,1,2,1,12)=0*import_infections;
            ydot6(2,1,1,1,2,1,1)=(0*import_infections)/max(1,(sum(sum(sum(sum(sum(sum(Y(2,:,:,:,:,:,12:20)))))))/max(1,sum(sum(sum(sum(sum(sum(Y(2,:,:,:,:,:,1:20)))))))))); % keep prevalence constant among former
            ydot6(3,1,1,1,2,1,12)=1*import_infections;
            
            if TT > 38 && TT < 66
                missing = interp1(dem(:,1), dem(:,2), TT+1950)... % population at time TT
                    - (sum(sum(sum(sum(sum(sum(Y(:,:,:,:,:,1,1:20),1),2),3),4),5),7) ... % minus total population
                    + sum(sum(sum(sum(sum(sum(ydot6(:,:,:,:,:,1,1:20),1),2),3),4),5),7)... % correct for other arrivals
                    + sum(sum(sum(sum(sum(sum(ydot3(:,:,:,:,:,1,1:20),1),2),3),4),5),7)); % correct for deaths
                ydot6(3,:,:,:,:,1,1) = ydot6(3,:,:,:,:,1,1) + ...
                    min(1/dt, missing ...
                    ./ sum(sum(sum(sum(sum(sum(Y(:,:,:,:,:,1,1)))))))).* Y(3,:,:,:,:,1,1);
            elseif TT >=66
                missing = interp1(dem(:,1), dem(:,2), 66+1950)... % population at time TT
                    - (sum(sum(sum(sum(sum(sum(Y(:,:,:,:,:,1,1:20),1),2),3),4),5),7) ... % minus total population
                    + sum(sum(sum(sum(sum(sum(ydot6(:,:,:,:,:,1,1:20),1),2),3),4),5),7)... % correct for other arrivals
                    + sum(sum(sum(sum(sum(sum(ydot3(:,:,:,:,:,1,1:20),1),2),3),4),5),7)); % correct for deaths
                ydot6(3,:,:,:,:,1,1) = ydot6(3,:,:,:,:,1,1) + ...
                    min(1/dt, missing ...
                    ./ sum(sum(sum(sum(sum(sum(Y(:,:,:,:,:,1,1)))))))).* Y(3,:,:,:,:,1,1);
            end
            
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
            progress = progression;
            
            if TT>=66 && TT <= 70 % scale-up of Antibody and RNA testing rates between now and 2020
                mult = min(1, 1 - (TT-66) / 4);
                progress(1,1,2,1) = mult*progression_base(1,1,2,1) + (1-mult)*progression(1,1,2,1);
                progress(2,1,2,1) = mult*progression_base(2,1,2,1) + (1-mult)*progression(2,1,2,1);
                progress(1,2,2,1) = mult*progression_base(1,2,2,1) + (1-mult)*progression(1,2,2,1);
                progress(2,2,2,1) = mult*progression_base(2,2,2,1) + (1-mult)*progression(2,2,2,1);
            end
            
            %only those in care progress due to zeroes in prog
            prog = permute(reshape(repmat(progress,1,num_age*num_intervention),num_pops,num_cascade,num_age,num_intervention,num_engagement,num_region),[1,2,3,4,5,6]);
            prog = min(1/dt, prog); % need to protect against adding too many people at once if time steps are too large
            
            prog(1:2,1,:,:,:,:) = (followup) * prog(1:2,1,:,:,:,:);
            %prog(1,:,:,:,1,:) = 0;
            if strcmp(scenario,'OST_test') == 1
                prog(1,1,:,1:2,:,:) = 0*prog(1,1,:,1:2,:,:); % In these scenarios don't test people who are not engaged in OST of NSP
            end
            ydot7=0*Y;
            
            ydot7(1,1:4,:,:,:,:,17:20) = -min(12, 1/dt)*Y(1,1:4,:,:,:,:,17:20); % assumes DC, HCC and LT HCV antibodies are detected and made ready for treatment immediately
            ydot7(1,6,:,:,2,:,17:20) = reshape(sum(sum(-ydot7(1,1:4,:,:,:,:,17:20),2),5),1,1,num_age,num_intervention, 1, num_region,4); % they all go to the second engagement level
            
            ydot7(:,[1,2,3,4],:,:,:,:,12) = ydot7(:,[1,2,3,4],:,:,:,:,12)-prog(:,[1,2,3,4],:,:,:,:).*Y(:,[1,2,3,4],:,:,:,:,12); % HCV antibody screening rate, RNA testing, genotype and liver testing for F0-F4 disease
            ydot7(:,[1,2,3,4],:,:,:,:,13) = ydot7(:,[1,2,3,4],:,:,:,:,13)-prog(:,[1,2,3,4],:,:,:,:).*Y(:,[1,2,3,4],:,:,:,:,13);
            ydot7(:,[1,2,3,4],:,:,:,:,14) = ydot7(:,[1,2,3,4],:,:,:,:,14)-prog(:,[1,2,3,4],:,:,:,:).*Y(:,[1,2,3,4],:,:,:,:,14);
            ydot7(:,[1,2,3,4],:,:,:,:,15) = ydot7(:,[1,2,3,4],:,:,:,:,15)-prog(:,[1,2,3,4],:,:,:,:).*Y(:,[1,2,3,4],:,:,:,:,15);
            ydot7(:,[1,2,3,4],:,:,:,:,16) = ydot7(:,[1,2,3,4],:,:,:,:,16)-prog(:,[1,2,3,4],:,:,:,:).*Y(:,[1,2,3,4],:,:,:,:,16);
            
            ydot7(:,[2,3,4,6],:,:,:,:,12) = ydot7(:,[2,3,4,6],:,:,:,:,12)+prog(:,[1,2,3,4],:,:,:,:).*Y(:,[1,2,3,4],:,:,:,:,12); % HCV antibody screening rate, RNA testing, genotype and liver testing for F0-F4 disease
            ydot7(:,[2,3,4,6],:,:,:,:,13) = ydot7(:,[2,3,4,6],:,:,:,:,13)+prog(:,[1,2,3,4],:,:,:,:).*Y(:,[1,2,3,4],:,:,:,:,13);
            ydot7(:,[2,3,4,6],:,:,:,:,14) = ydot7(:,[2,3,4,6],:,:,:,:,14)+prog(:,[1,2,3,4],:,:,:,:).*Y(:,[1,2,3,4],:,:,:,:,14);
            ydot7(:,[2,3,4,6],:,:,:,:,15) = ydot7(:,[2,3,4,6],:,:,:,:,15)+prog(:,[1,2,3,4],:,:,:,:).*Y(:,[1,2,3,4],:,:,:,:,15);
            ydot7(:,[2,3,4,6],:,:,:,:,16) = ydot7(:,[2,3,4,6],:,:,:,:,16)+prog(:,[1,2,3,4],:,:,:,:).*Y(:,[1,2,3,4],:,:,:,:,16);
            
            ydot7(:,2,:,:,:,:,28) = +prog(:,1,:,:,:,:).*sum(Y(:,1,:,:,:,:,12:16),7) + sum(min(12,1/dt)*Y(:,1,:,:,:,:,17:20),7);%Total antibody tests. Assumed to occur at the same rate for all PWID
            ydot7(:,3,:,:,:,:,29) = +prog(:,2,:,:,:,:).*sum(Y(:,2,:,:,:,:,12:16),7) + sum(sum(min(12,1/dt)*Y(:,1:2,:,:,:,:,17:20),2),7); % Total RNA detections
            ydot7(:,4,:,:,:,:,30) = +prog(:,3,:,:,:,:).*sum(Y(:,3,:,:,:,:,12:16),7) + sum(sum(min(12,1/dt)*Y(:,1:3,:,:,:,:,17:20),2),7); % total genotypes completed
            ydot7(:,6,:,:,:,:,31) = +prog(:,4,:,:,:,:).*sum(Y(:,4,:,:,:,:,12:15),7); % total liver tests (fibroscans) early disease stage
            ydot7(:,6,:,:,:,:,32) = +prog(:,4,:,:,:,:).*Y(:,4,:,:,:,:,16)...
                + reshape(sum(sum(min(12,1/dt)*Y(:,1:4,:,:,:,:,17:20),2),7),num_pops,1,num_age,num_intervention, num_engagement, num_region); % total liver tests (fibroscans) late disease stage
            

            if TT > 67 && APRI == 1 % From 2017 fibroscan not required for F2 or less if scenario uses APRI. Also genotype removed
                if cascade_scale_time==0 mult=1; else mult = min(1, 1 - ((cascade_scale_time - (TT-67)) / cascade_scale_time)); end
                prog2 = prog;
                prog2(:,4,:,:,2:end,:) = ((1-mult)*prog2(:,4,:,:,2:end,:)+mult*52.*prog2(:,4,:,:,2:end,:)./prog2(:,4,:,:,2:end,:)); % cascade progression rates super fast meaning skipped fibroscan
                prog2 = min(1/dt, prog2);
                ydot7(:,4,:,:,2:end,:,12) = prog2(:,3,:,:,2:end,:).*Y(:,3,:,:,2:end,:,12)- prog2(:,4,:,:,2:end,:).*Y(:,4,:,:,2:end,:,12);  % Fibroscan skipped for F0-F2
                ydot7(:,4,:,:,2:end,:,13) = prog2(:,3,:,:,2:end,:).*Y(:,3,:,:,2:end,:,13)- prog2(:,4,:,:,2:end,:).*Y(:,4,:,:,2:end,:,13);
                ydot7(:,4,:,:,2:end,:,14) = prog2(:,3,:,:,2:end,:).*Y(:,3,:,:,2:end,:,14)- prog2(:,4,:,:,2:end,:).*Y(:,4,:,:,2:end,:,14);
                
                ydot7(:,6,:,:,2:end,:,12) = prog2(:,4,:,:,2:end,:).*Y(:,4,:,:,2:end,:,12); % ready for treatment
                ydot7(:,6,:,:,2:end,:,13) = prog2(:,4,:,:,2:end,:).*Y(:,4,:,:,2:end,:,13);
                ydot7(:,6,:,:,2:end,:,14) = prog2(:,4,:,:,2:end,:).*Y(:,4,:,:,2:end,:,14);
                
                ydot7(:,6,:,:,:,:,31) = +prog(:,4,:,:,:,:).*Y(:,4,:,:,:,:,15);
                
            end
            
            
            
            %% Treatment
            ydot8=0*Y;
            % f has zeroes for people not in care.
            ydot8(:,6,:,:,:,:,12) = - min(f(:,6,:,:,:,:,12), prog(:,5,:,:,:,:).*Y(:,6,:,:,:,:,12)); % treatment from F0 - LT2
            ydot8(:,6,:,:,:,:,13) = - min(f(:,6,:,:,:,:,13), prog(:,5,:,:,:,:).*Y(:,6,:,:,:,:,13));
            ydot8(:,6,:,:,:,:,14) = - min(f(:,6,:,:,:,:,14), prog(:,5,:,:,:,:).*Y(:,6,:,:,:,:,14));
            ydot8(:,6,:,:,:,:,15) = - min(f(:,6,:,:,:,:,15), prog(:,5,:,:,:,:).*Y(:,6,:,:,:,:,15));
            ydot8(:,6,:,:,:,:,16) = - min(f(:,6,:,:,:,:,16), prog(:,5,:,:,:,:).*Y(:,6,:,:,:,:,16));
            ydot8(:,6,:,:,:,:,17) = - min(f(:,6,:,:,:,:,17), prog(:,5,:,:,:,:).*Y(:,6,:,:,:,:,17));
            ydot8(:,6,:,:,:,:,18) = - min(f(:,6,:,:,:,:,18), prog(:,5,:,:,:,:).*Y(:,6,:,:,:,:,18));
            ydot8(:,6,:,:,:,:,19) = - min(f(:,6,:,:,:,:,19), prog(:,5,:,:,:,:).*Y(:,6,:,:,:,:,19));
            ydot8(:,6,:,:,:,:,20) = - min(f(:,6,:,:,:,:,20), prog(:,5,:,:,:,:).*Y(:,6,:,:,:,:,20));
            
            completion = rate(alpha*p_complete,alpha);
            
            ydot8(:,7,:,:,:,:,7) = completion(:,6,:,:,:,:).*(-ydot8(:,6,:,:,:,:,12)); % in treatment from F0 and going to succeed
            ydot8(:,7,:,:,:,:,8) = completion(:,6,:,:,:,:).*(-ydot8(:,6,:,:,:,:,13)); % T1
            ydot8(:,7,:,:,:,:,9) = completion(:,6,:,:,:,:).*(-ydot8(:,6,:,:,:,:,14)); % T2
            ydot8(:,7,:,:,:,:,10) = completion(:,6,:,:,:,:).*(-ydot8(:,6,:,:,:,:,15)); % T3
            ydot8(:,7,:,:,:,:,11) = completion(:,6,:,:,:,:).*sum((-ydot8(:,6,:,:,:,:,16:20)),7); % T4
            
            ydot8(:,8,:,:,:,:,12) = (1-completion(:,6,:,:,:,:)).*(-ydot8(:,6,:,:,:,:,12)); % failing treatment
            ydot8(:,8,:,:,:,:,13) = (1-completion(:,6,:,:,:,:)).*(-ydot8(:,6,:,:,:,:,13));
            ydot8(:,8,:,:,:,:,14) = (1-completion(:,6,:,:,:,:)).*(-ydot8(:,6,:,:,:,:,14));
            ydot8(:,8,:,:,:,:,15) = (1-completion(:,6,:,:,:,:)).*(-ydot8(:,6,:,:,:,:,15));
            ydot8(:,8,:,:,:,:,16) = (1-completion(:,6,:,:,:,:)).*(-ydot8(:,6,:,:,:,:,16));
            ydot8(:,8,:,:,:,:,17) = (1-completion(:,6,:,:,:,:)).*(-ydot8(:,6,:,:,:,:,17));
            ydot8(:,8,:,:,:,:,18) = (1-completion(:,6,:,:,:,:)).*(-ydot8(:,6,:,:,:,:,18));
            ydot8(:,8,:,:,:,:,19) = (1-completion(:,6,:,:,:,:)).*(-ydot8(:,6,:,:,:,:,19));
            ydot8(:,8,:,:,:,:,20) = (1-completion(:,6,:,:,:,:)).*(-ydot8(:,6,:,:,:,:,20));
            
            ydot8(:,7,:,:,2:end,:,7:11) = ydot8(:,7,:,:,2:end,:,7:11) - min(omega, 1/dt) *Y(:,7,:,:,2:end,:,7:11); % Recovering and moving to the 1-st level. Must be engaged in care
            ydot8(:,1,:,:,2:end,:,1:5) = + min(omega, 1/dt)*Y(:,7,:,:,2:end,:,7:11); % Recovering and moving to the 1-st level
            
            ydot8(:,6,:,:,:,:,23) = - sum(ydot8(:,6,:,:,:,:,12:15),7); % total first round treatments from F0-F3
            ydot8(:,6,:,:,:,:,25) = - sum(ydot8(:,6,:,:,:,:,16:20),7); % total first round treatments from F4 onwards
            
            
            %% Retreatment
            ydot9=0*Y;
            ydot9(:,8,:,:,:,:,12) = - min(f(:,8,:,:,:,:,12), prog(:,6,:,:,:,:).*Y(:,8,:,:,:,:,12)); % retreatment from F0 - LT2
            ydot9(:,8,:,:,:,:,13) = - min(f(:,8,:,:,:,:,13), prog(:,6,:,:,:,:).*Y(:,8,:,:,:,:,13));
            ydot9(:,8,:,:,:,:,14) = - min(f(:,8,:,:,:,:,14), prog(:,6,:,:,:,:).*Y(:,8,:,:,:,:,14));
            ydot9(:,8,:,:,:,:,15) = - min(f(:,8,:,:,:,:,15), prog(:,6,:,:,:,:).*Y(:,8,:,:,:,:,15));
            ydot9(:,8,:,:,:,:,16) = - min(f(:,8,:,:,:,:,16), prog(:,6,:,:,:,:).*Y(:,8,:,:,:,:,16));
            ydot9(:,8,:,:,:,:,17) = - min(f(:,8,:,:,:,:,17), prog(:,6,:,:,:,:).*Y(:,8,:,:,:,:,17));
            ydot9(:,8,:,:,:,:,18) = - min(f(:,8,:,:,:,:,18), prog(:,6,:,:,:,:).*Y(:,8,:,:,:,:,18));
            ydot9(:,8,:,:,:,:,19) = - min(f(:,8,:,:,:,:,19), prog(:,6,:,:,:,:).*Y(:,8,:,:,:,:,19));
            ydot9(:,8,:,:,:,:,20) = - min(f(:,8,:,:,:,:,20), prog(:,6,:,:,:,:).*Y(:,8,:,:,:,:,20));
            
            ydot9(:,9,:,:,:,:,7) = completion(:,8,:,:,:,:).*(-ydot9(:,8,:,:,:,:,12)); % in treatment from F0 and going to succeed
            ydot9(:,9,:,:,:,:,8) = completion(:,8,:,:,:,:).*(-ydot9(:,8,:,:,:,:,13)); % T1
            ydot9(:,9,:,:,:,:,9) = completion(:,8,:,:,:,:).*(-ydot9(:,8,:,:,:,:,14)); % T2
            ydot9(:,9,:,:,:,:,10) = completion(:,8,:,:,:,:).*(-ydot9(:,8,:,:,:,:,15)); % T3
            ydot9(:,9,:,:,:,:,11) = completion(:,8,:,:,:,:).*sum((-ydot9(:,8,:,:,:,:,16:20)),7); % T4
            
            ydot9(:,10,:,:,:,:,12) = (1-completion(:,8,:,:,:,:)).*(-ydot9(:,8,:,:,:,:,12)); % failing treatment
            ydot9(:,10,:,:,:,:,13) = (1-completion(:,8,:,:,:,:)).*(-ydot9(:,8,:,:,:,:,13));
            ydot9(:,10,:,:,:,:,14) = (1-completion(:,8,:,:,:,:)).*(-ydot9(:,8,:,:,:,:,14));
            ydot9(:,10,:,:,:,:,15) = (1-completion(:,8,:,:,:,:)).*(-ydot9(:,8,:,:,:,:,15));
            ydot9(:,10,:,:,:,:,16) = (1-completion(:,8,:,:,:,:)).*(-ydot9(:,8,:,:,:,:,16));
            ydot9(:,10,:,:,:,:,17) = (1-completion(:,8,:,:,:,:)).*(-ydot9(:,8,:,:,:,:,17));
            ydot9(:,10,:,:,:,:,18) = (1-completion(:,8,:,:,:,:)).*(-ydot9(:,8,:,:,:,:,18));
            ydot9(:,10,:,:,:,:,19) = (1-completion(:,8,:,:,:,:)).*(-ydot9(:,8,:,:,:,:,19));
            ydot9(:,10,:,:,:,:,20) = (1-completion(:,8,:,:,:,:)).*(-ydot9(:,8,:,:,:,:,20));
            
            ydot9(:,9,:,:,2:end,:,7:11) = ydot9(:,9,:,:,2:end,:,7:11) - min(omega, 1/dt)*Y(:,9,:,:,2:end,:,7:11); % Recovering and moving to the 1-st level. Must be engaged in care
            ydot9(:,1,:,:,2:end,:,1:5) = + min(omega, 1/dt)*Y(:,9,:,:,2:end,:,7:11); % Recovering and moving to the 1-st level
            
            ydot9(:,8,:,:,:,:,23) = - sum(ydot9(:,8,:,:,:,:,12:15),7); % total first round treatments from F0-F3
            ydot9(:,8,:,:,:,:,25) = - sum(ydot9(:,8,:,:,:,:,16:20),7); % total first round treatments from F4 onwards
            
            %% Engagement in care
            ydot10 = 0*Y;
            ydot10(:,:,:,:,1,:,16:20) = -min(12, 1/dt)*Y(:,:,:,:,1,:,16:20); % sick people engaged in care
            ydot10(:,:,:,:,2,:,16:20) = min(12, 1/dt)*Y(:,:,:,:,1,:,16:20);
            
            ydot10(:,:,:,:,1,:,:) = ydot10(:,:,:,:,1,:,:) - min(1*exit_IDU, 1/dt)*Y(:,:,:,:,1,:,:)...
                + min(1*exit_IDU, 1/dt)*Y(:,:,:,:,2,:,:); 
            ydot10(:,:,:,:,2,:,:) = ydot10(:,:,:,:,2,:,:) + min(1*exit_IDU, 1/dt)*Y(:,:,:,:,1,:,:)...
                - min(1*exit_IDU, 1/dt)*Y(:,:,:,:,2,:,:);
            

            %% Changes in intervention coverage
            ydot11 = 0*Y;
            
            
            ydot11(1,:,:,1,:,:,1:20) = ydot11(1,:,:,1,:,:,1:20) + ...
                min(1/dt, ((1-nsp_coverage)*(1-ost_coverage)*sum(Y(1,:,:,:,:,:,1:20),4) - Y(1,:,:,1,:,:,1:20))./sum(Y(1,:,:,:,:,:,1:20),4))...
                .*Y(1,:,:,1,:,:,1:20);
            ydot11(1,:,:,2,:,:,1:20) = ydot11(1,:,:,2,:,:,1:20) + ...
                min(1/dt, (nsp_coverage*(1-ost_coverage)*sum(Y(1,:,:,:,:,:,1:20),4) - Y(1,:,:,2,:,:,1:20))./sum(Y(1,:,:,:,:,:,1:20),4))...
                .*Y(1,:,:,2,:,:,1:20);
            ydot11(1,:,:,3,:,:,1:20) = ydot11(1,:,:,3,:,:,1:20) + ...
                min(1/dt, ((1-nsp_coverage)*ost_coverage*sum(Y(1,:,:,:,:,:,1:20),4) - Y(1,:,:,3,:,:,1:20))./sum(Y(1,:,:,:,:,:,1:20),4))...
                .*Y(1,:,:,3,:,:,1:20);
            ydot11(1,:,:,4,:,:,1:20) = ydot11(1,:,:,4,:,:,1:20) + ...
                min(1/dt, (nsp_coverage*ost_coverage*sum(Y(1,:,:,:,:,:,1:20),4) - Y(1,:,:,4,:,:,1:20))./sum(Y(1,:,:,:,:,:,1:20),4))...
                .*Y(1,:,:,4,:,:,1:20);

                
            %OST
            ydot11(1,:,:,1:2,:,:,1:20) = ydot11(1,:,:,1:2,:,:,1:20) - min(1/ost_duration, 1/dt) * Y(1,:,:,1:2,:,:,1:20);
            ydot11(1,:,:,3:4,:,:,1:20) = ydot11(1,:,:,3:4,:,:,1:20) + min(1/ost_duration, 1/dt) * Y(1,:,:,1:2,:,:,1:20);
            ydot11(1,:,:,3:4,:,:,1:20) = ydot11(1,:,:,3:4,:,:,1:20) - min(1/ost_duration, 1/dt) * Y(1,:,:,3:4,:,:,1:20);
            ydot11(1,:,:,1:2,:,:,1:20) = ydot11(1,:,:,1:2,:,:,1:20) + min(1/ost_duration, 1/dt) * Y(1,:,:,3:4,:,:,1:20);
            
%             ydot11(2:3,:,:,3:4,:,:,1:20) = ydot11(2:3,:,:,3:4,:,:,1:20) - min(1/ost_duration, 1/dt) * Y(2:3,:,:,3:4,:,:,1:20); %former and other can drop out of OST but not recruit
%             ydot11(2:3,:,:,1:2,:,:,1:20) = ydot11(2:3,:,:,1:2,:,:,1:20) + min(1/ost_duration, 1/dt) * Y(2:3,:,:,3:4,:,:,1:20);
            
            %NSP
            ydot11(1,:,:,[1,3],:,:,1:20) = ydot11(1,:,:,[1,3],:,:,1:20) - min(1/nsp_duration, 1/dt) * Y(1,:,:,[1,3],:,:,1:20);
            ydot11(1,:,:,[2,4],:,:,1:20) = ydot11(1,:,:,[2,4],:,:,1:20) + min(1/nsp_duration, 1/dt) * Y(1,:,:,[1,3],:,:,1:20);
            ydot11(1,:,:,[2,4],:,:,1:20) = ydot11(1,:,:,[2,4],:,:,1:20) - min(1/nsp_duration, 1/dt) * Y(1,:,:,[2,4],:,:,1:20);
            ydot11(1,:,:,[1,3],:,:,1:20) = ydot11(1,:,:,[1,3],:,:,1:20) + min(1/nsp_duration, 1/dt) * Y(1,:,:,[2,4],:,:,1:20);
%             
%             if TT > 66 && harm_reduction_coverage >= 0
%                 ydot11(1,:,:,1,:,:,1:20) = min(1/dt,sum((1-harm_reduction_coverage)*Y(1,:,:,:,:,:,1:20),4) - Y(1,:,:,1,:,:,1:20));
%                 ydot11(1,:,:,2,:,:,1:20) = min(1/dt,1/3*sum(harm_reduction_coverage*Y(1,:,:,:,:,:,1:20),4) - Y(1,:,:,2,:,:,1:20));
%                 ydot11(1,:,:,3,:,:,1:20) = min(1/dt,1/3*sum(harm_reduction_coverage*Y(1,:,:,:,:,:,1:20),4) - Y(1,:,:,3,:,:,1:20));
%                 ydot11(1,:,:,4,:,:,1:20) = min(1/dt,1/3*sum(harm_reduction_coverage*Y(1,:,:,:,:,:,1:20),4) - Y(1,:,:,4,:,:,1:20));
%             end
            

            %% Combine
            ydot=real(reshape(ydot1+ydot3+ydot4+ydot5+ydot6+ydot7+ydot8+ydot9+ydot10+ydot11,num_pops*num_cascade*num_age*num_intervention*num_engagement*num_region*(27+6),1));
            %clear Y ydot1 ydot3 ydot4 ydot5 ydot6 ydot7 ydot8 ydot9 ydot10 ydot11 prog MU lambda2 ...
            %    S S1 S2 S3 S4 A T T1 T2 T3 T4 F0 F1 F2 F3 F4 DC HCC LT LT2 F4_transfer Ldeath T_total HCC_transfer T_F4on_total Liver_transplants Inc Cas1 Cas2 Cas3 Cas4 Cas5 Cas6
            
            
        else
            ydot = zeros(num_pops*num_cascade*num_age*num_intervention*num_engagement*num_region*(27+6),1);
        end
        %% Rate function
        
        
        y(t,:)=y(t-1,:)+max(dt*ydot', -y(t-1,:));
    end
    function ad = rate(r_XY_PWID,r_XY) % sub-function that distributes factors for PWID and everyone else
        ad = reshape(repmat([r_XY_PWID;r_XY;r_XY],num_cascade,num_age,num_intervention,num_engagement,num_region),num_pops,num_cascade,num_age,num_intervention, num_engagement,num_region);
    end
y=reshape(y,size(y,1),num_pops,num_cascade,num_age,num_intervention,num_engagement,num_region,(27+6));
TT=tvec;%TT(1:end);
TT=[t0(1:end-1);TT];

y1=y;
for i=2:length(y) %Make sure these are the additional treatments/HCC transfers/LTs/liver transplants/incidence, rather than total
    y1(i,:,:,:,:,:,:,[21,23:(27+6)])=max(y(i,:,:,:,:,:,:,[21,23:(27+6)])-y(i-1,:,:,:,:,:,:,[21,23:(27+6)]),0);
end
y=y1;


end

