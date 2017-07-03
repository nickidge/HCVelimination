function [R0] = R0(t,y_in)

global mu_PWID mu_former exit_IDU r_relapse delta alpha p_complete omega infect total_PWID PWID0...
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2...
    c1_chronicPWID c2_PWID c3_PWID c4_PWID c5_PWID c6_PWID ...
    c1_chronicformerPWID c2_formerPWID c3_formerPWID c4_formerPWID c5_formerPWID c6_formerPWID ...
    c1_chronic c2 c3 c4 c5 c6 imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9 ...
    scenario cascade_scale_time age_mix start_year treat ...
    APRI num_pops num_cascade num_age num_intervention num_engagement num_region infect_factor progression


Fty = Sim_sub2(t,y_in);
thresh = 10^(-6) * ones(num_pops*num_cascade*num_age*num_intervention*num_engagement*num_region*20,1);
fac = [];
vectorized = 1;
[dFdy,fac] = numjac('Sim_sub2',t,y_in,Fty,thresh,fac,vectorized);%,S,g)

Fty2 = Sim_sub3(t,y_in);
fac2 = [];
[dFdy2,fac2] = numjac('Sim_sub3',t,y_in,Fty2,thresh,fac,vectorized);%,S,g)

R0 = max(abs(eig(dFdy2 * inv(dFdy))));

    function[R0_out]=Sim_sub3(TT,y)
        Y=reshape(y,num_pops, num_cascade, num_age, num_intervention, num_engagement,num_region,20);
        
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

        pop=reshape(sum(sum(sum(sum(sum(Y(:,:,:,:,:,:,1:20),2),3),4),5),7),num_pops, num_region); %vector of populaiton sizes
        
        R_inc=1;
        
        %Version that gets ~9% in F34 at 2015
        if TT<40 R_inc=1+ (1.5/40)*(TT); elseif TT>=40 && TT<45 R_inc=2.5-(1.5/5)*(TT-40); elseif TT>=45 R_inc=1; end
        
        function ad = rate(r_XY_PWID,r_XY) % sub-function that distributes factors for PWID and everyone else
            ad = reshape(repmat([r_XY_PWID;r_XY;r_XY;r_XY;r_XY_PWID;r_XY],num_cascade,num_age,num_intervention,num_engagement,num_region),num_pops,num_cascade,num_age,num_intervention, num_engagement,num_region);
        end
        
        % Force of infection
        if age_mix < 0
            lambda = zeros(num_pops,num_region);
            for i = 1:num_pops
                for j = 1:num_region
                    if pop(i,j) > 0
                        if i==1 || i==5
                            p = [1,5];
                        elseif i==4 || i==5 || i==6
                            p = [4,5,6];
                        else
                            p = i;
                        end
                        lambda(i,j)=R_inc*infect*sum(sum(sum(sum(sum(A(p,:,:,:,:,j)+F0(p,:,:,:,:,j)+F1(p,:,:,:,:,j)+F2(p,:,:,:,:,j)+F3(p,:,:,:,:,j)+F4(p,:,:,:,:,j)+DC(p,:,:,:,:,j)+HCC(p,:,:,:,:,j)+LT(p,:,:,:,:,j)+LT2(p ,:,:,:,:,j))))))/sum(pop(p,j));
                    else
                        lambda(i,j)=0;
                    end
                end
            end
            
            lambda2 =  permute(reshape(repmat(lambda,1,num_age),num_pops,num_region,num_age),[1,3,2]);
        else
            lambda_age = zeros(num_pops,num_age,num_region);
            for j = 1:num_region
                for i =1:num_pops
                    if i==1 % identifying which population groups mix
                        p = [1,5]; % all PWID population groups
                    elseif i==4 || i==6
                        p = [4,5,6]; % all MSM population groups
                    elseif i ==5 % MSM and PWID
                        p = [1,4,5,6];
                    else
                        p = i;
                    end
                    for a=2:num_age-1
                        if sum(sum(sum(sum(sum(Y(i,:,a,:,:,:,1:20)))))) > 0
                            lambda_age(i,a,j) = max(0, R_inc*infect*(...
                                age_mix * sum(sum(sum(sum(A(p,:,a,:,:,j)+F0(p,:,a,:,:,j)+F1(p,:,a,:,:,j)+F2(p,:,a,:,:,j)+F3(p,:,a,:,:,j)+F4(p,:,a,:,:,j)+DC(p,:,a,:,:,j)+HCC(p,:,a,:,:,j)+LT(p,:,a,:,:,j)+LT2(p,:,a,:,:,j)))))/sum(sum(sum(sum(sum(Y(p,:,a,:,:,j,1:20)))))) + ...
                                (1-age_mix)*sum(sum(sum(sum(sum(A(p,:,[1:a-1,a+1:end],:,:,j)+F0(p,:,[1:a-1,a+1:end],:,:,j)+F1(p,:,[1:a-1,a+1:end],:,:,j)+F2(p,:,[1:a-1,a+1:end],:,:,j)+F3(p,:,[1:a-1,a+1:end],:,:,j)+F4(p,:,[1:a-1,a+1:end],:,:,j)+DC(p,:,[1:a-1,a+1:end],:,:,j)+HCC(p,:,[1:a-1,a+1:end],:,:,j)+LT(p,:,[1:a-1,a+1:end],:,:,j)+LT2(p,:,[1:a-1,a+1:end],:,:,j))))))/sum(sum(sum(sum(sum(sum(Y(p,:,[1:a-1,a+1:end],:,:,j,1:20))))))) ...
                                ));
                        end
                    end
                    if sum(sum(sum(sum(Y(i,:,1,:,:,j,1:20))))) > 0 && sum(sum(sum(sum(sum(Y(i,:,2:end,:,:,j,1:20)))))) > 0
                        lambda_age(i,1,j) = max(0,R_inc*infect*(...
                            age_mix*sum(sum(sum(sum(A(p,:,1,:,:,j)+F0(p,:,1,:,:,j)+F1(p,:,1,:,:,j)+F2(p,:,1,:,:,j)+F3(p,:,1,:,:,j)+F4(p,:,1,:,:,j)+DC(p,:,1,:,:,j)+HCC(p,:,1,:,:,j)+LT(p,:,1,:,:,j)+LT2(p,:,1,:,:,j)))))/sum(sum(sum(sum(sum(Y(p,:,1,:,:,j,1:20)))))) + ...
                            (1-age_mix)*sum(sum(sum(sum(sum(A(p,:,2:end,:,:,j)+F0(p,:,2:end,:,:,j)+F1(p,:,2:end,:,:,j)+F2(p,:,2:end,:,:,j)+F3(p,:,2:end,:,:,j)+F4(p,:,2:end,:,:,j)+DC(p,:,2:end,:,:,j)+HCC(p,:,2:end,:,:,j)+LT(p,:,2:end,:,:,j)+LT2(p,:,2:end,:,:,j))))))/sum(sum(sum(sum(sum(sum(Y(p,:,2:end,:,:,j,1:20))))))) ...
                            ));
                    end
                    if sum(sum(sum(Y(i,:,num_age,1:20)))) > 0
                        lambda_age(i,num_age,j) = max(0, R_inc*infect*(...
                            age_mix*sum(sum(sum(sum(A(p,:,num_age,:,:,j)+F0(p,:,num_age,:,:,j)+F1(p,:,num_age,:,:,j)+F2(p,:,num_age,:,:,j)+F3(p,:,num_age,:,:,j)+F4(p,:,num_age,:,:,j)+DC(p,:,num_age,:,:,j)+HCC(p,:,num_age,:,:,j)+LT(p,:,num_age,:,:,j)+LT2(p,:,num_age,:,:,j)))))/sum(sum(sum(sum(sum(Y(p,:,num_age,:,:,j,1:20)))))) + ...
                            (1-age_mix)*sum(sum(sum(sum(sum(A(p,:,1:end-1,:,:,j)+F0(p,:,1:end-1,:,:,j)+F1(p,:,1:end-1,:,:,j)+F2(p,:,1:end-1,:,:,j)+F3(p,:,1:end-1,:,:,j)+F4(p,:,1:end-1,:,:,j)+DC(p,:,1:end-1,:,:,j)+HCC(p,:,1:end-1,:,:,j)+LT(p,:,1:end-1,:,:,j)+LT2(p,:,1:end-1,:,:,j))))))/sum(sum(sum(sum(sum(sum(Y(p,:,1:end-1,:,:,j,1:20))))))) ...
                            ));
                    end
                end
            end
            lambda_age = permute(reshape(repmat(infect_factor,1,num_age*num_cascade*num_engagement),num_pops,num_intervention, num_region,num_age, num_cascade,num_engagement),[1,5,4,2,6,3]) ...
                .* permute(reshape(repmat(lambda_age,1,num_cascade*num_engagement*num_intervention),num_pops,num_age,num_region,num_cascade,num_intervention,num_engagement),[1,4,2,5,6,3]);
            lambda2 = lambda_age;
        end
        
        R0_a = zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement,num_region,20);
        R0_a(1,1,:,6) = lambda2(1,1,:,:,:,:)*(1-delta).*S(1,1,:,:,:,:);
        R0_a(1,1,:,13) = lambda2(1,1,:,:,:,:).*S1(1,1,:,:,:,:);
        R0_a(1,1,:,14) = lambda2(1,1,:,:,:,:).*S2(1,1,:,:,:,:);
        R0_a(1,1,:,15) = lambda2(1,1,:,:,:,:).*S3(1,1,:,:,:,:);
        R0_a(1,1,:,16) = lambda2(1,1,:,:,:,:).*S4(1,1,:,:,:,:);
        
        R0_out = reshape(R0_a, num_pops*num_cascade*num_age*num_intervention*num_engagement*num_region*20,1);
    end

%ODE subfunction
    function[ydot]=Sim_sub2(TT,y)
        
        Y=reshape(y,num_pops, num_cascade, num_age, num_intervention, num_engagement,num_region,20);
        
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
        
        pop=reshape(sum(sum(sum(sum(sum(Y(:,:,:,:,:,:,1:20),2),3),4),5),7),num_pops, num_region); %vector of populaiton sizes
        
        R_inc=1;
        
        %Version that gets ~9% in F34 at 2015
        if TT<40 R_inc=1+ (1.5/40)*(TT); elseif TT>=40 && TT<45 R_inc=2.5-(1.5/5)*(TT-40); elseif TT>=45 R_inc=1; end
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
        
        function ad = rate(r_XY_PWID,r_XY) % sub-function that distributes factors for PWID and everyone else
            ad = reshape(repmat([r_XY_PWID;r_XY;r_XY;r_XY;r_XY_PWID;r_XY],num_cascade,num_age,num_intervention,num_engagement,num_region),num_pops,num_cascade,num_age,num_intervention, num_engagement,num_region);
        end
        
        if TT > start_year
            % Force of infection
            if age_mix < 0
                lambda = zeros(num_pops,num_region);
                for i = 1:num_pops
                    for j = 1:num_region
                        if pop(i,j) > 0
                            if i==1 || i==5
                                p = [1,5];
                            elseif i==4 || i==5 || i==6
                                p = [4,5,6];
                            else
                                p = i;
                            end
                            lambda(i,j)=R_inc*infect*sum(sum(sum(sum(sum(A(p,:,:,:,:,j)+F0(p,:,:,:,:,j)+F1(p,:,:,:,:,j)+F2(p,:,:,:,:,j)+F3(p,:,:,:,:,j)+F4(p,:,:,:,:,j)+DC(p,:,:,:,:,j)+HCC(p,:,:,:,:,j)+LT(p,:,:,:,:,j)+LT2(p ,:,:,:,:,j))))))/sum(pop(p,j));
                        else
                            lambda(i,j)=0;
                        end
                    end
                end
                
                lambda2 =  permute(reshape(repmat(lambda,1,num_age),num_pops,num_region,num_age),[1,3,2]);
            else
                lambda_age = zeros(num_pops,num_age,num_region);
                for j = 1:num_region
                    for i =1:num_pops
                        if i==1 % identifying which population groups mix
                            p = [1,5]; % all PWID population groups
                        elseif i==4 || i==6
                            p = [4,5,6]; % all MSM population groups
                        elseif i ==5 % MSM and PWID
                            p = [1,4,5,6];
                        else
                            p = i;
                        end
                        for a=2:num_age-1
                            if sum(sum(sum(sum(sum(Y(i,:,a,:,:,:,1:20)))))) > 0
                                lambda_age(i,a,j) = max(0, R_inc*infect*(...
                                    age_mix * sum(sum(sum(sum(A(p,:,a,:,:,j)+F0(p,:,a,:,:,j)+F1(p,:,a,:,:,j)+F2(p,:,a,:,:,j)+F3(p,:,a,:,:,j)+F4(p,:,a,:,:,j)+DC(p,:,a,:,:,j)+HCC(p,:,a,:,:,j)+LT(p,:,a,:,:,j)+LT2(p,:,a,:,:,j)))))/sum(sum(sum(sum(sum(Y(p,:,a,:,:,j,1:20)))))) + ...
                                    (1-age_mix)*sum(sum(sum(sum(sum(A(p,:,[1:a-1,a+1:end],:,:,j)+F0(p,:,[1:a-1,a+1:end],:,:,j)+F1(p,:,[1:a-1,a+1:end],:,:,j)+F2(p,:,[1:a-1,a+1:end],:,:,j)+F3(p,:,[1:a-1,a+1:end],:,:,j)+F4(p,:,[1:a-1,a+1:end],:,:,j)+DC(p,:,[1:a-1,a+1:end],:,:,j)+HCC(p,:,[1:a-1,a+1:end],:,:,j)+LT(p,:,[1:a-1,a+1:end],:,:,j)+LT2(p,:,[1:a-1,a+1:end],:,:,j))))))/sum(sum(sum(sum(sum(sum(Y(p,:,[1:a-1,a+1:end],:,:,j,1:20))))))) ...
                                    ));
                            end
                        end
                        if sum(sum(sum(sum(Y(i,:,1,:,:,j,1:20))))) > 0 && sum(sum(sum(sum(sum(Y(i,:,2:end,:,:,j,1:20)))))) > 0
                            lambda_age(i,1,j) = max(0,R_inc*infect*(...
                                age_mix*sum(sum(sum(sum(A(p,:,1,:,:,j)+F0(p,:,1,:,:,j)+F1(p,:,1,:,:,j)+F2(p,:,1,:,:,j)+F3(p,:,1,:,:,j)+F4(p,:,1,:,:,j)+DC(p,:,1,:,:,j)+HCC(p,:,1,:,:,j)+LT(p,:,1,:,:,j)+LT2(p,:,1,:,:,j)))))/sum(sum(sum(sum(sum(Y(p,:,1,:,:,j,1:20)))))) + ...
                                (1-age_mix)*sum(sum(sum(sum(sum(A(p,:,2:end,:,:,j)+F0(p,:,2:end,:,:,j)+F1(p,:,2:end,:,:,j)+F2(p,:,2:end,:,:,j)+F3(p,:,2:end,:,:,j)+F4(p,:,2:end,:,:,j)+DC(p,:,2:end,:,:,j)+HCC(p,:,2:end,:,:,j)+LT(p,:,2:end,:,:,j)+LT2(p,:,2:end,:,:,j))))))/sum(sum(sum(sum(sum(sum(Y(p,:,2:end,:,:,j,1:20))))))) ...
                                ));
                        end
                        if sum(sum(sum(Y(i,:,num_age,1:20)))) > 0
                            lambda_age(i,num_age,j) = max(0, R_inc*infect*(...
                                age_mix*sum(sum(sum(sum(A(p,:,num_age,:,:,j)+F0(p,:,num_age,:,:,j)+F1(p,:,num_age,:,:,j)+F2(p,:,num_age,:,:,j)+F3(p,:,num_age,:,:,j)+F4(p,:,num_age,:,:,j)+DC(p,:,num_age,:,:,j)+HCC(p,:,num_age,:,:,j)+LT(p,:,num_age,:,:,j)+LT2(p,:,num_age,:,:,j)))))/sum(sum(sum(sum(sum(Y(p,:,num_age,:,:,j,1:20)))))) + ...
                                (1-age_mix)*sum(sum(sum(sum(sum(A(p,:,1:end-1,:,:,j)+F0(p,:,1:end-1,:,:,j)+F1(p,:,1:end-1,:,:,j)+F2(p,:,1:end-1,:,:,j)+F3(p,:,1:end-1,:,:,j)+F4(p,:,1:end-1,:,:,j)+DC(p,:,1:end-1,:,:,j)+HCC(p,:,1:end-1,:,:,j)+LT(p,:,1:end-1,:,:,j)+LT2(p,:,1:end-1,:,:,j))))))/sum(sum(sum(sum(sum(sum(Y(p,:,1:end-1,:,:,j,1:20))))))) ...
                                ));
                        end
                    end
                end
                lambda_age = permute(reshape(repmat(infect_factor,1,num_age*num_cascade*num_engagement),num_pops,num_intervention, num_region,num_age, num_cascade,num_engagement),[1,5,4,2,6,3]) ...
                    .* permute(reshape(repmat(lambda_age,1,num_cascade*num_engagement*num_intervention),num_pops,num_age,num_region,num_cascade,num_intervention,num_engagement),[1,4,2,5,6,3]);
                lambda2 = lambda_age;
            end
            
            %TREATMENT ALLOCATION FUNCTION, f
            
            f=100*ones(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region,20);
            if alpha < 0.8 && TT > 52
                c5_PWID = -log(1-0.02*sum(sum(sum(sum(sum(sum(Y(1,[6,8],:,:,2:end,:,12:20)))))))/sum(sum(sum(sum(sum(sum(sum(Y(1,:,:,:,2:end,:,1:20))))))))); % tear 2% of infected p.a.
                c5_formerPWID = -log(1-0.02*sum(sum(sum(sum(sum(sum(sum(Y(:,[6,8],:,:,2:end,:,12:20))))))))/sum(sum(sum(sum(sum(sum(sum(Y(:,:,:,:,2:end,:,1:20)))))))));
                c5 = -log(1-0.02*sum(sum(sum(sum(sum(sum(sum(Y(:,[6,8],:,:,2:end,:,12:20))))))))/sum(sum(sum(sum(sum(sum(sum(Y(:,:,:,:,2:end,:,1:20)))))))));
                c6_PWID = 0; % no retreatment with IFN
                c6_formerPWID = 0;
                c6 = 0;
            end
            if TT>52 && TT<=66
                trea=treat(3);
                trea_extra=0;
                [f]=treatment_alloc(Y,trea,trea_extra,0);
            end
            if TT>66
                trea=treat(1);
                if TT<=67 trea_extra=treat(2); else trea_extra=0; end
                [f]=treatment_alloc(Y,trea,trea_extra,0);
            end
            
            if TT>=66 imported2 = rate(imported,0); else imported2 = rate(0,0); end
            new0 = (lambda2+imported2).*S;
            new1 = (lambda2+imported2).*S1;
            new2 = (lambda2+imported2).*S2;
            new3 = (lambda2+imported2).*S3;
            new4 = (lambda2+imported2).*S4;
            
            
            %Infection and progression
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
                r_LT1LT2*LT-r_LTdeath2*LT2); %LT2
            
            
            %Mortality
            MU = permute(reshape(repmat([mu_PWID;mu_former;mu_former;mu_former;mu_PWID;mu_former],1,num_cascade*num_intervention*num_engagement*num_region),...
                num_pops,num_age,num_cascade,num_intervention,num_engagement,num_region),[1,3,2,4,5,6]); % assumes current and former die at same rate
            l_death=sum(sum(sum(sum(sum(sum(r_HCCdeath*HCC+r_LTdeath1*LT+r_LTdeath2*LT2+r_DCdeath*DC+r_S4death*S4+r_S4death*T4))))));
            ydot3 = 0*Y;
            ydot3(:,:,:,:,:,:,1:20)=cat(7,-S.*MU,...
                -S1.*MU,-S2.*MU,-S3.*MU,-S4.*MU,-A.*MU,-T.*MU,-T1.*MU,-T2.*MU,-T3.*MU,-T4.*MU,...
                -F0.*MU,-F1.*MU,-F2.*MU,-F3.*MU,-F4.*MU,-DC.*MU,-HCC.*MU,-LT.*MU,-LT2.*MU);
            
            
            %Cessation and relapse into IDU
            r_relapse_rel=r_relapse;%*R_inc;
            exit_IDU_rel=exit_IDU;%/R_inc;
            
            ydot4=0*Y;
            ydot4([1,5],:,:,:,:,:,1:20)=-Y([1,5],:,:,:,:,:,1:20)*exit_IDU_rel + Y([2,6],:,:,:,:,:,1:20)*r_relapse;
            ydot4([2,6],:,:,:,:,:,1:20)=Y([1,5],:,:,:,:,:,1:20)*exit_IDU_rel-Y([2,6],:,:,:,:,:,1:20)*r_relapse- 1/(10*12)*Y([2,6],:,:,:,:,:,1:20); % former for 10 years
            ydot4([3,4],:,:,:,:,:,1:20)=1/(10*12)*Y([2,6],:,:,:,:,:,1:20);
            
            
            %Let people age
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
            
            %New arrivals
            IDUs=PWID0;
            if TT<50
                IDUs = (TT-start_year)*(100000-PWID0)/(50-start_year) + PWID0;
            elseif TT>=50 && TT<55
                IDUs = 100000 - (TT-50)*(100000 -total_PWID)/5;
            elseif TT>=55
                IDUs = total_PWID;
            end
            mu_arrivals = max(IDUs-sum(sum(sum(sum(sum(sum(sum(Y([1,5],:,:,:,:,:,1:20))))))))...
                +sum(sum(sum(sum(sum(sum(sum(-ydot3([1,5],:,:,:,:,:,:))))))))+sum(sum(sum(sum(sum(sum(sum(-ydot4([1,5],:,:,:,:,:,1:20)))))))),...
                0);
            %         mu_arrivals=sum(sum(sum(sum(-ydot3(:,:,:,:)))));
            ydot6=0*Y;
            ydot6(1,1,1,1,1,1,1)=l_death+mu_arrivals;
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
            elseif TT>=65 import_infections=imp9; end
            ydot6(2,1,1,1,1,1,12)=0.5*import_infections;
            ydot6(2,1,1,1,1,1,1)=(0.5*import_infections)/max(1,(sum(sum(sum(sum(sum(sum(Y(2,:,:,:,:,:,12:20)))))))/max(1,sum(sum(sum(sum(sum(sum(Y(2,:,:,:,:,:,1:20)))))))))); % keep prevalence constant among former
            ydot6(3,1,1,1,1,1,12)=0.5*import_infections;
            
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
            prog = permute(reshape(repmat(progression,1,num_age*num_intervention),num_pops,num_cascade,num_engagement,num_region,num_age,num_intervention),[1,2,5,6,3,4]);
            %only those in care progress due to zeroes in prog
            ydot7=0*Y;
            
            ydot7(1,1:4,:,:,:,:,17:20) = -12*Y(1,1:4,:,:,:,:,17:20); % assumes DC, HCC and LT HCV antibodies are detected and made ready for treatment immediately
            ydot7(1,6,:,:,2,:,17:20) = reshape(sum(sum(12*Y(1,1:4,:,:,:,:,17:20),2),5),1,1,num_age,num_intervention, 1, num_region,4); % they all go to the second engagement level
            
            ydot7(:,[1,2,3,4],:,:,:,:,12:16) = -reshape(repmat(prog(:,[1,2,3,4],:,:,:,:),1,5), num_pops, 4, num_age, num_intervention, num_engagement, num_region,5)...
                .*Y(:,[1,2,3,4],:,:,:,:,12:16); % HCV antibody screening rate, RNA testing, genotype and liver testing for F0-F4 disease
            ydot7(:,[2,3,4,6],:,:,:,:,12:16) = +reshape(repmat(prog(:,[1,2,3,4],:,:,:,:),1,5), num_pops, 4, num_age, num_intervention, num_engagement, num_region,5)...
                .*Y(:,[2,3,4,6],:,:,:,:,12:16); % HCV antibody screening rate, RNA testing, genotype and liver testing for F0-F4 disease
            
            
            if TT > 67 && APRI == 1 % Fibroscan not required for F2 or less if scenario uses APRI
                if cascade_scale_time==0 mult=1; else mult = min(1, 1 - ((cascade_scale_time - (TT-67)) / cascade_scale_time)); end
                prog2 = prog;
                prog2(:,4,:,:,2:end,:) = ((1-mult)*prog2(:,4,:,:,2:end,:)+mult*52);
                ydot7(:,4,:,:,2:end,:,12:14) = reshape(repmat(prog2(:,3,:,:,2:end,:),1,3),num_pops,1,num_age, num_intervention,num_engagement-1,num_region,3)...
                    .*Y(:,3,:,:,2:end,:,12:14)...
                    - reshape(repmat(prog2(:,4,:,:,2:end,:),1,3),num_pops,1,num_age, num_intervention,num_engagement-1,num_region,3)...
                    .*Y(:,4,:,:,2:end,:,12:14); % Fibroscan skipped for F0-F2
                ydot7(:,6,:,:,2:end,:,12:14) = reshape(repmat(prog2(:,4,:,:,2:end,:),1,3),num_pops,1,num_age, num_intervention,num_engagement-1,num_region,3)...
                    .*Y(:,4,:,:,2:end,:,12:14); % ready for treatment

            end
            
            %Treatment
            ydot8=0*Y;
            % f has zeroes for people not in care.
            ydot8(:,6,:,:,:,:,12:20) = - min(f(:,6,:,:,:,:,12:20), ...
                reshape(repmat(prog(:,5,:,:,:,:),1,9), num_pops, 1, num_age, num_intervention, num_engagement, num_region,9)...
                .*Y(:,6,:,:,:,:,12:20)); % treatment from F0 - LT2
            
            completion = rate(alpha*p_complete,alpha);
            
            ydot8(:,7,:,:,:,:,7) = completion(:,6,:,:,:,:).*(-ydot8(:,6,:,:,:,:,12)); % in treatment from F0 and going to succeed
            ydot8(:,7,:,:,:,:,8) = completion(:,6,:,:,:,:).*(-ydot8(:,6,:,:,:,:,13)); % T1
            ydot8(:,7,:,:,:,:,9) = completion(:,6,:,:,:,:).*(-ydot8(:,6,:,:,:,:,14)); % T2
            ydot8(:,7,:,:,:,:,10) = completion(:,6,:,:,:,:).*(-ydot8(:,6,:,:,:,:,15)); % T3
            ydot8(:,7,:,:,:,:,11) = completion(:,6,:,:,:,:).*sum((-ydot8(:,6,:,:,:,:,16:20)),7); % T4
            
            ydot8(:,8,:,:,:,:,12:20) = reshape(repmat(1-completion(:,8,:,:,:,:),1,9),num_pops, 1, num_age, num_intervention, num_engagement, num_region,9)...
                .*(-ydot8(:,6,:,:,:,:,12:20)); % failing treatment
            
            ydot8(:,7,:,:,2:end,:,7:11) = ydot8(:,7,:,:,2:end,:,7:11) - omega*Y(:,7,:,:,2:end,:,7:11); % Recovering and moving to the 1-st level. Must be engaged in care
            ydot8(:,1,:,:,2:end,:,7:11) = + omega*Y(:,7,:,:,2:end,:,7:11); % Recovering and moving to the 1-st level

            
            
            %Retreatment
            ydot9=0*Y;
            ydot9(:,8,:,:,:,:,12:20) = - min(f(:,8,:,:,:,:,12:20), ...
                reshape(repmat(prog(:,5,:,:,:,:),1,9), num_pops, 1, num_age, num_intervention, num_engagement, num_region,9)...
                .*Y(:,8,:,:,:,:,12:20)); % treatment from F0 - LT2
            
            ydot9(:,9,:,:,:,:,7) = completion(:,8,:,:,:,:).*(-ydot9(:,8,:,:,:,:,12)); % in treatment from F0 and going to succeed
            ydot9(:,9,:,:,:,:,8) = completion(:,8,:,:,:,:).*(-ydot9(:,8,:,:,:,:,13)); % T1
            ydot9(:,9,:,:,:,:,9) = completion(:,8,:,:,:,:).*(-ydot9(:,8,:,:,:,:,14)); % T2
            ydot9(:,9,:,:,:,:,10) = completion(:,8,:,:,:,:).*(-ydot9(:,8,:,:,:,:,15)); % T3
            ydot9(:,9,:,:,:,:,11) = completion(:,8,:,:,:,:).*sum((-ydot9(:,8,:,:,:,:,16:20)),7); % T4
            
            ydot9(:,10,:,:,:,:,12:20) = reshape(repmat(1-completion(:,8,:,:,:,:),1,9),num_pops, 1, num_age, num_intervention, num_engagement, num_region,9)...
                .*(-ydot9(:,8,:,:,:,:,12:20)); % failing treatment
            
            ydot9(:,9,:,:,2:end,:,7:11) = ydot9(:,9,:,:,2:end,:,7:11) - omega*Y(:,9,:,:,2:end,:,7:11); % Recovering and moving to the 1-st level. Must be engaged in care
            ydot9(:,1,:,:,2:end,:,7:11) = + omega*Y(:,9,:,:,2:end,:,7:11); % Recovering and moving to the 1-st level
            
            
            % Engagement in care
            ydot10 = 0*Y;
            
            
            % Incarceration
            ydot11 = 0*Y;
            
            ydot=real(reshape(ydot1+ydot3+ydot4+ydot5+ydot6+ydot7+ydot8+ydot9,num_pops*num_cascade*num_age*num_intervention*num_engagement*num_region*(20),1));
            
            
        else
            ydot = zeros(num_pops*num_cascade*num_age*num_intervention*num_engagement*num_region*(20),1);
        end
        
    end
end