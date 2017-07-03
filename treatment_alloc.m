function[f]=treatment_alloc(Y, trea, trea_extra, cal)
global target_late num_pops num_cascade num_age num_intervention num_engagement num_region
%% Treatment
t=[6,8]; % treatments only occur in t-th level of cascade. Treatments also only go to those engaged in care
if cal == 1 
    target_late=1;
    num_pops = 3;
    num_cascade = 10;
    num_age = 9;
    num_intervention = 1;
    num_engagement = 2;
    num_region = 1;
elseif cal==0
    if target_late<0 target_late=max(0,sum(sum(sum(sum(sum(sum(sum(Y(:,t,:,:,2:end,:,15:20))))))))/sum(sum(sum(sum(sum(sum(sum(Y(:,t,:,:,2:end,:,12:20))))))))); end %Proportion of infected who are in late liver disease stage (i.e to treat proportionally)
end
tot_infected_PWID=reshape(sum(sum(sum(sum(sum(sum(Y(1,t,:,:,2:end,:,12:20))))))),1,1);
tot_PWID_adv=reshape(sum(sum(sum(sum(sum(sum(Y(1,t,:,:,2:end,:,15:20))))))),1,1);
tot_nonPWID_adv=reshape(sum(sum(sum(sum(sum(sum(sum(Y(2:end,t,:,:,2:end,:,15:20)))))))),1,1);
%[tot_infected_PWID,tot_PWID_adv, tot_nonPWID_adv, target_late]
intersection = max(0,tot_PWID_adv/(target_late*(tot_PWID_adv + tot_nonPWID_adv) + (1-target_late)*tot_infected_PWID)); %The proportion of treatments to be allocated to infected PWID with late liver disease

f=zeros(num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region,20);
trea = trea+trea_extra;
if intersection*trea <= tot_PWID_adv && tot_PWID_adv > 0 %PWID with advanced liver disease
    f(1,t,:,:,2:end,:,15:20)=intersection*trea*Y(1,t,:,:,2:end,:,15:20)/tot_PWID_adv;
else
    f(1,t,:,:,2:end,:,15:20)=Y(1,t,:,:,2:end,:,15:20);
end
if target_late*(1-intersection)*trea <= tot_nonPWID_adv && tot_nonPWID_adv > 0 %Former PWID with advanced liver disease
    f(2:end,t,:,:,2:end,:,15:20)=target_late*(1-intersection)*trea*Y(2:end,t,:,:,2:end,:,15:20)/tot_nonPWID_adv;
else
    f(2:end,t,:,:,2:end,:,15:20)=Y(2:end,t,:,:,2:end,:,15:20);
end
if (1-target_late)*(1-intersection)*trea <= tot_infected_PWID-tot_PWID_adv && tot_infected_PWID-tot_PWID_adv > 0 %PWID with early liver disease
    f(1,t,:,:,2:end,:,12:14)=(1-target_late)*(1-intersection)*trea*Y(1,t,:,:,2:end,:,12:14)/(tot_infected_PWID-tot_PWID_adv);
else
    f(1,t,:,:,2:end,:,12:14)=Y(1,t,:,:,2:end,:,12:14);
end
denom = sum(sum(sum(sum(sum(sum(sum(Y(:,t,:,:,2:end,:,12:20)-f(:,t,:,:,2:end,:,12:20)))))))); % leftover people to treat
treat_remaining = trea - sum(sum(sum(sum(sum(sum(sum(f(:,t,:,:,2:end,:,12:20))))))));
if treat_remaining>0 && treat_remaining<=denom
    f(:,t,:,:,2:end,:,12:20)=f(:,t,:,:,2:end,:,12:20)+treat_remaining*(Y(:,t,:,:,2:end,:,12:20)-f(:,t,:,:,2:end,:,12:20))/denom;
elseif treat_remaining>0 && treat_remaining>denom
    f(:,t,:,:,2:end,:,12:20)=Y(:,t,:,:,2:end,:,12:20);
end
f=max(0,f);
end

