function[R0_out]=Sim_sub3(TT,y)
global mu_PWID mu_former exit_IDU r_relapse delta alpha p_complete omega infect total_PWID PWID0...
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2...
    c1_chronicPWID c2_PWID c3_PWID c4_PWID c5_PWID c6_PWID ...
    c1_chronicformerPWID c2_formerPWID c3_formerPWID c4_formerPWID c5_formerPWID c6_formerPWID ...
    c1_chronic c2 c3 c4 c5 c6 imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9 ...
    scenario cascade_scale_time age_mix start_year ...
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



R0_a = zeros(3,10,9,20);
R0_a(1,1,:,6) = lambda2(1,1,:)*(1-delta).*S(1,1,:);
R0_a(1,1,:,13) = lambda2(1,1,:).*S1(1,1,:);
R0_a(1,1,:,14) = lambda2(1,1,:).*S2(1,1,:);
R0_a(1,1,:,15) = lambda2(1,1,:).*S3(1,1,:);
R0_a(1,1,:,16) = lambda2(1,1,:).*S4(1,1,:);

R0_out = reshape(R0_a, 3*10*9*20,1);
end
