function [ycomb_noage, summary, tr, tr_] = gather_outputs(y1,y2,TT)
global Tin num_pops num_cascade num_age num_intervention num_engagement num_region

ycomb=[y1(1:end-1,:,:,:,:,:,:,:);y2];
%Make everyone susceptible for R0 calculation
ycomb_S = zeros(length(TT),num_pops, num_cascade, num_age, num_intervention, num_engagement, num_region,20);
ycomb_S(:,:,1,:,:,:,:,1) = ycomb(:,:,1,:,:,:,:,1) + sum(sum(ycomb(:,:,2:10,:,:,:,:,[6,7,12]),3),8);
ycomb_S(:,:,1,:,:,:,:,2) = ycomb(:,:,1,:,:,:,:,2) + sum(sum(ycomb(:,:,2:10,:,:,:,:,[8,13]),3),8);
ycomb_S(:,:,1,:,:,:,:,3) = ycomb(:,:,1,:,:,:,:,3) + sum(sum(ycomb(:,:,2:10,:,:,:,:,[9,14]),3),8);
ycomb_S(:,:,1,:,:,:,:,4) = ycomb(:,:,1,:,:,:,:,4) + sum(sum(ycomb(:,:,2:10,:,:,:,:,[10,15]),3),8);
ycomb_S(:,:,1,:,:,:,:,5) = ycomb(:,:,1,:,:,:,:,5) + sum(sum(ycomb(:,:,2:10,:,:,:,:,[11,16:20]),3),8);

y2_noage=reshape(sum(sum(sum(sum(y2,4),5),6),7),size(y2,1),num_pops, num_cascade,27+6); %Reshape to sum over age stratification, interventions, engagement and region
ycomb_noage=reshape(sum(sum(sum(sum(ycomb,4),5),6),7),size(ycomb,1),num_pops, num_cascade,27+6); %Reshape to sum over age stratification
%y3=squeeze(sum(sum(sum(y2,4),6),7)); %Reshape to sum over age stratification, interventions, engagement and region
Tint=Tin+[0,1,13]; t_val=zeros(1,length(Tint));
for l=1:length(Tint) %Find times corresponding to years since intervention
    t_val(l)=find(TT>=(Tint(l)),1);
end
TT_=TT(t_val(2):t_val(3))-Tint(2);
[cost,DALY,life_exp,tr]=Costs_age(y2_noage([t_val(2):t_val(3)]-t_val(1)+1,:,:,:),TT_,y2([t_val(2):t_val(3)]-t_val(1)+1,:,:,:,:,:,:,:));
tr_(1) = sum(sum(sum(y2_noage([t_val(2):t_val(3)]-t_val(1)+1,1,:,[23,25])))); % Total treatments for PWID
tr_(2) = sum(sum(sum(y2_noage([t_val(2):t_val(3)]-t_val(1)+1,2,:,[23,25])))); % Total treatments for former PWID
tr_(3) = sum(sum(sum(y2_noage([t_val(2):t_val(3)]-t_val(1)+1,3,:,[23,25])))); % Total treatments for other
tr_(4) = sum(sum(sum(sum(y2_noage([t_val(2):t_val(3)]-t_val(1)+1,:,:,[23,25]))))); % Total treatments

%Put a bunch of summary statistics together
summary = [cost,DALY,tr(4),tr_(4),...
    sum(sum(sum(ycomb_noage(find(TT>=80,1)+1:find(TT>=81,1),:,:,27))))/(TT(find(TT>=81,1))-TT(find(TT>=80,1))),... %Incidence in 2030
    (sum(sum(ycomb_noage(find(TT>=81,1),1:3,:,22)))-sum(sum(ycomb_noage(find(TT>=80,1),1:3,:,22))))...
    /(TT(find(TT>=81,1))-TT(find(TT>=80,1))),... %Liver related deaths in 2030
    100*sum(sum(ycomb_noage(t_val(1),1,:,[6,12:20])))./sum(sum(ycomb_noage(t_val(1),1,:,1:20))),... %Prevalence among PWID in 2015
    100*sum(sum(ycomb_noage(t_val(3),1,:,[6,12:20])))./sum(sum(ycomb_noage(t_val(3),1,:,1:20))),... %Prevalence among PWID in 2030
    sum(sum(ycomb_noage(t_val(3),1:3,:,22))),... %Liver related deaths 2015 to 2030
    100*sum(sum(ycomb_noage(find(TT>=69,1),1,:,[6,12:20])))./sum(sum(ycomb_noage(find(TT>=69,1),1,:,1:20))),... %Prevalence among PWID in 2019
    100*sum(sum(ycomb_noage(find(TT>=70,1),1,:,[6,12:20])))./sum(sum(ycomb_noage(find(TT>=70,1),1,:,1:20))),... %Prevalence among PWID in 2020
    100*sum(sum(ycomb_noage(find(TT>=75,1),1,:,[6,12:20])))./sum(sum(ycomb_noage(find(TT>=75,1),1,:,1:20)))]; %Prevalence among PWID in 2025

end