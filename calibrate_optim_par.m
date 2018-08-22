function [output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases,output_ost, output_nsp, output_diagnoses] = calibrate_optim(maxiter, swarmsize)
global infect prev0 progression imp1 imp2 imp3 imp4 imp5 imp6 imp7 imp8 imp9 cascade0 cascade0_PWID disease0 cases0 ost0 nsp0 diagnoses0...
    data y0_init t0_init y0 t0 treat Tin infect_base progression_base ost_enrollment nsp_enrollment r_inc_up start_year followup...
    num_pops num_cascade num_age num_intervention num_engagement num_region infect_factor treat_projected ...
    r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID ost_coverage nsp_coverage 

data = {prev0, cascade0, cascade0_PWID, disease0, cases0, ost0, nsp0, diagnoses0}; % prevalence, cascade to calibrate to
lag = 1; % iterations between successive samples
xp = [10*infect,...
    progression(1,1:4,2,1), ... %PWID first 4 steps
    progression(2,1:4,2,1),... %former PWID first 4 steps
    progression(3,1:4,2,1),... % Other first 4 steps
    imp1, imp2, imp3, imp4, imp5, imp6, imp7, imp8, imp9, ...
    ost_enrollment, nsp_enrollment,r_inc_up, start_year,...
    r_F0F1, r_F1F2, r_F2F3, r_F3F4, r_F0F1_PWID, r_F1F2_PWID, r_F2F3_PWID, r_F3F4_PWID];

nvars = length(xp);

lb = max(0.1*xp, 0.001*ones(1,length(xp)));
lb(1) = 0.01; lb(14:16) = 1000; lb(22) = 0; 
lb(23) = 0; lb(24) = 0; 
lb(25) = 0; lb(26) = 0; % lower bounds for height of rel_incidence function and epidemic start year
lb(27:34) = 0.5*xp(27:34);
ub = 10*xp;
ub(14:19) = 10000; ub(20:22) = 2000; %ub(22) = 0; 
ub(23) = 0.1; ub(24) = 0.1; % Forcing no NSP or OST
ub(25) = 5; ub(26) = 40; % upper bounds for height of rel_incidence function and epidemic start year
ub(27:34) = 1.5*xp(27:34);

%options.UseVectorized = true;
ms = MultiStart('UseParallel', true);
options = optimoptions('particleswarm', 'Display', 'iter', 'MaxIter', maxiter,'StallIterlimit',20,'Swarmsize',swarmsize,'UseParallel',false);
[x, fval, exitflag, output] = particleswarm(@objective,nvars,lb,ub, options)

[infect, progression, imp1, imp2, imp3, imp4, imp5 ,imp6 ,imp7, imp8, imp9,...
    ost_enrollment, nsp_enrollment, r_inc_up, start_year,...
    r_F0F1, r_F1F2, r_F2F3, r_F3F4, r_F0F1_PWID, r_F1F2_PWID, r_F2F3_PWID, r_F3F4_PWID] = feval(@assign,x);

Run=16; %Years to run the model

P=1000; %Population to start the model
infected0=0.10; % initial proportion of PWID infected (in 1950)
PWID0=P*0.5; %Equilibrium proportion of PWID to former PWID
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

for k = 6:13
    if x(k) < 0.9*x(k-4)
        x(k) = x(k-4);
    end
end
for k = 17:22
    if x(k) > 1.1*x(k-1)
        x(k) = 1.1*x(k-1);
    end
end
% for k = 15:(length(x)-2)
%     if x(k)< 0.9*x(k-1)
%         x(k) = x(k-1);
%     end
% end

infect = x(1) / 10;
progression(1,1:4,2,1) = x(2:5);
progression(2,1:4,2,1) = x(6:9);
progression(3,1:4,2,1) = x(10:13);
imp1 = x(14);
imp2 = x(15);
imp3 = x(16);
imp4 = x(17);
imp5 = x(18);
imp6 = x(19);
imp7 = x(20);
imp8 = x(21);
imp9 = x(22);
ost_enrollment = x(23);
nsp_enrollment = x(24);
r_inc_up = x(25);
start_year = x(26);
r_F0F1 = x(27);
r_F1F2 = x(28);
r_F2F3 = x(29);
r_F3F4  = x(30);
r_F0F1_PWID  = x(31);
r_F1F2_PWID  = x(32);
r_F2F3_PWID  = x(33);
r_F3F4_PWID = x(34);

nsp_coverage = nsp0(1,2);
ost_coverage = ost0(1,2);

[TT, y] = DE_track_age(Tin, y0, 0, treat);

[output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases, output_ost, output_nsp, output_diagnoses] = model_vals(TT,y)



%%

    function probX = objective(xp)
        
        %loaddata
        
        [infect, progression, imp1, imp2, imp3, imp4, imp5 ,imp6 ,imp7, imp8, imp9, ost_enrollment, nsp_enrollment, r_inc_up, start_year,...
            r_F0F1, r_F1F2, r_F2F3, r_F3F4, r_F0F1_PWID, r_F1F2_PWID, r_F2F3_PWID, r_F3F4_PWID] = feval(@assign,xp);
        
         
        data = {prev0, cascade0, cascade0_PWID, disease0, cases0, ost0, nsp0, diagnoses0}; % prevalence, cascade to calibrate to
        
        nsp_coverage = nsp0(1,2);
        ost_coverage = ost0(1,2);
        
        for k = 6:13
            if xp(k) < 0.9*xp(k-4)
                xp(k) = xp(k-4);
            end
        end
        for k = 17:22
            if xp(k) > 1.1*xp(k-1)
                xp(k) = 1.1*xp(k-1);
            end
        end
%         for k = 15:(length(xp)-2)
%             if xp(k)< 0.9*xp(k-1)
%                 xp(k) = xp(k-1);
%             end
%         end
        
        infect = xp(1) / 10;
        progression(1,1:4,2,1) = xp(2:5);
        progression(2,1:4,2,1) = xp(6:9);
        progression(3,1:4,2,1) = xp(10:13);
        imp1 = xp(14);
        imp2 = xp(15);
        imp3 = xp(16);
        imp4 = xp(17);
        imp5 = xp(18);
        imp6 = xp(19);
        imp7 = xp(20);
        imp8 = xp(21);
        imp9 = xp(22);
        ost_enrollment = xp(23);
        nsp_enrollment = xp(24);
        r_inc_up = xp(25);
        start_year = xp(26);
        [imp1, imp2, imp3, imp4, imp5, imp6, imp7, imp8 imp9];
        r_F0F1 = xp(27);
        r_F1F2 = xp(28);
        r_F2F3 = xp(29);
        r_F3F4  = xp(30);
        r_F0F1_PWID  = xp(31);
        r_F1F2_PWID  = xp(32);
        r_F2F3_PWID  = xp(33);
        r_F3F4_PWID = xp(34);
                
        probX = 0;
        
        
        [TT, y] = DE_track_age(Tin, y0_init, 0, treat);
        %infect, progression, imp1, imp2, imp3, imp4, imp5 ,imp6 ,imp7);
        
        [output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases, output_ost, output_nsp, output_diagnoses] = model_vals(TT,y);

        sigma_prev = 0.001*data{1}(:,2:end); %standard deviations for prevalence by year
        sigma_cascade = 0.1*data{2}(:,2:end); %standard deviations for cascade by year
        sigma_cascade_PWID = 0.1*data{3}(:,2:end); %standard deviations for cascade by year
        sigma_disease = 0.01*data{4}(:,2:end); %standard deviations for disease by year
        sigma_cases = 0.001*data{5}(:,2:end); %standard deviations for cases by year
        sigma_ost = 0.1*data{6}(:,2:end); %standard deviations for proportion of PWID on ost by year
        sigma_nsp = 0.1*data{7}(:,2:end); %standard deviations for proportion of PWID accessing NSP by year
        sigma_diagnoses = 1*data{8}(:,2:end); %standard deviations for proportion of PWID accessing NSP by year
        
        for years = 1:length(output_prev(:,1))
            probX = probX - (-(data{1}(years,2)-output_prev(years))^2/(2*sigma_prev(years)^2)); %*...
            %exp(-(data(2)-output(2))^2/(2*sigma(2)^2))/sigma(2); % Gaussian with SD
        end
%         for years = 1:length(output_cascade(:,1))
%             for vals = 1:length(output_cascade(years,:))
%                 probX = probX - (-(data{2}(years,vals+1)-output_cascade(years, vals))^2/(2*sigma_cascade(vals)^2));%/sigma_cascade(vals); %*...
%             end
%         end
%         for years = 1:length(output_cascade_PWID(:,1))
%             for vals = 1:length(output_cascade_PWID(years,:))
%                 probX = probX - (-(data{3}(years,vals+1)-output_cascade_PWID(years, vals))^2/(2*sigma_cascade_PWID(vals)^2));%/sigma_cascade(vals); %*...
%             end
%         end
        if data{4}(1,1) ~= 0
            for years = 1:length(output_disease(:,1))
                for vals = 1:length(output_disease(years,:))
                    probX = probX - (-(data{4}(years,vals+1)-output_disease(years,vals))^2/(2*sigma_disease(years)^2));%/sigma_disease(years); %*...
                end
            end
        end
        for years = 1:length(output_cases(:,1))
            probX = probX - (-(data{5}(years,2)-output_cases(years))^2/(2*sigma_cases(years)^2));%/sigma_cases(years); %*...
        end
%         for years = 1:length(output_ost(:,1))
%             probX = probX - (-(data{6}(years,2)-output_ost(years))^2/(2*sigma_ost(years)^2));%/sigma_cases(years); %*...
%         end
%         for years = 1:length(output_nsp(:,1))
%             probX = probX - (-(data{7}(years,2)-output_nsp(years))^2/(2*sigma_nsp(years)^2));%/sigma_cases(years); %*...
%         end
%         for years = 1:length(output_diagnoses(:,1))
%             probX = probX - (-(data{8}(years,2)-output_diagnoses(years))^2/(2*sigma_diagnoses(years)^2));%/sigma_cases(years); %*...
%         end
    end

    function [output_prev, output_cascade, output_cascade_PWID, output_disease, output_cases, output_ost, output_nsp, output_diagnoses] = model_vals(TT,y)
        output_prev = zeros(length(prev0(:,1)),1);
        output_cascade = zeros(length(cascade0(:,1)),1);
        output_cascade_PWID = zeros(length(cascade0_PWID(:,1)),1);
        output_disease = zeros(length(disease0(:,1)),1);
        output_cases = zeros(length(cases0(:,1)),1);
        output_ost = zeros(length(ost0(:,1)),1);
        output_nsp = zeros(length(nsp0(:,1)),1);
        output_diagnoses = zeros(length(diagnoses0(:,1)),1);
        
        for years = 1:length(prev0(:,1))
            output_prev(years) = sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-prev0(years,1)),1),1,:,:,:,:,1,6:20))))))./sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-prev0(years,1)),1),1,:,:,:,:,1,1:20))))));
        end
        for years = 1:length(cascade0(:,1))
            output_cascade(years,1) = sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),:,2:10,:,:,:,1,6:20)))))))./sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),:,1:10,:,:,:,1,6:20)))))));
            output_cascade(years,2) = sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),:,3:10,:,:,:,1,6:20)))))))./sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),:,1:10,:,:,:,1,6:20)))))));
            output_cascade(years,3) = sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),:,4:10,:,:,:,1,6:20)))))))./sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),:,1:10,:,:,:,1,6:20)))))));
            output_cascade(years,4) = sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),:,6:10,:,:,:,1,6:20)))))))./sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),:,1:10,:,:,:,1,6:20)))))));
            %output_cascade(years,5) = sum(sum(sum(y(find(TT>=Tin-(2015-cascade0(years,1)),1),:,8,:,6:20))))./sum(sum(sum(sum(y(find(TT>=Tin-(2015-cascade0(years,1)),1),:,:,:,6:20)))));
        end
        for years = 1:length(cascade0_PWID(:,1))
            output_cascade_PWID(years,1) = sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0_PWID(years,1)),1),1,2:10,:,:,:,1,6:20))))))./sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),1,1:10,:,:,:,1,6:20))))));
            output_cascade_PWID(years,2) = sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0_PWID(years,1)),1),1,3:10,:,:,:,1,6:20))))))./sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),1,1:10,:,:,:,1,6:20))))));
            output_cascade_PWID(years,3) = sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0_PWID(years,1)),1),1,4:10,:,:,:,1,6:20))))))./sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),1,1:10,:,:,:,1,6:20))))));
            output_cascade_PWID(years,4) = sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0_PWID(years,1)),1),1,6:10,:,:,:,1,6:20))))))./sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cascade0(years,1)),1),1,1:10,:,:,:,1,6:20))))));
            %output_cascade(years,5) = sum(sum(sum(y(find(TT>=Tin-(2015-cascade0(years,1)),1),:,8,:,6:20))))./sum(sum(sum(sum(y(find(TT>=Tin-(2015-cascade0(years,1)),1),:,:,:,6:20)))));
        end
        
        for years = 1:length(disease0(:,1))
            dummy = max(2016, disease0(years,1));
            output_disease(years,1) = sum(sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-dummy),1),:,:,:,:,:,:,12:13))))))))./sum(sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-dummy),1),:,:,:,:,:,:,12:20)))))))); %F0-1
            output_disease(years,2) = sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-dummy),1),:,:,:,:,:,:,14)))))))./sum(sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-dummy),1),:,:,:,:,:,:,12:20)))))))); %F2
            output_disease(years,3) = sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-dummy),1),:,:,:,:,:,:,15)))))))./sum(sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-dummy),1),:,:,:,:,:,:,12:20)))))))); %F3
            output_disease(years,4) = sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-dummy),1),:,:,:,:,:,:,16)))))))./sum(sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-dummy),1),:,:,:,:,:,:,12:20)))))))); %F4
            output_disease(years,5) = sum(sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-dummy),1),:,:,:,:,:,:,17:20))))))))./sum(sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-dummy),1),:,:,:,:,:,:,12:20))))))));
        end
        for years = 1:length(cases0(:,1))
            output_cases(years) = sum(sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-cases0(years,1)),1),:,1:10,:,:,:,:,[6:20]))))))));
        end
        for years = 1:length(ost0(:,1))
            output_ost(years) = sum(sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-ost0(years,1)),1),1,:,:,3:4,:,:,1:20)))))))) / sum(sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-ost0(years,1)),1),1,:,:,:,:,:,1:20))))))));
        end
        for years = 1:length(nsp0(:,1))
            output_nsp(years) = sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-nsp0(years,1)),1),1,:,:,[2,4],:,1,1:20))))))) / sum(sum(sum(sum(sum(sum(y(find(TT>=Tin-(2016-nsp0(years,1)),1),1,:,:,:,:,1,1:20)))))));
        end
        for years = 1:length(diagnoses0(:,1))
            temp = Tin-(2016-diagnoses0(years,1));
            output_diagnoses(years) = sum(sum(sum(sum(sum(sum(sum(...
                y(find(TT>=temp-1,1):find(TT>=temp,1)-1,:,2,:,:,:,:,28))))))));
        end
    end

    function[infect, progression, imp1, imp2, imp3, imp4, imp5 ,imp6 ,imp7, imp8, imp9, ost_enrollment, nsp_enrollment, r_inc_up, start_year,...
            r_F0F1, r_F1F2, r_F2F3, r_F3F4, r_F0F1_PWID, r_F1F2_PWID, r_F2F3_PWID, r_F3F4_PWID] = assign(x)
        %loaddata
        infect = x(1) / 10;
        progression = zeros(num_pops,num_cascade,num_engagement,num_region);
        progression(1,1:4,2,1) = x(2:5);
        progression(2,1:4,2,1) = x(6:9);
        progression(3,1:4,2,1) = x(10:13);

        imp1 = x(14);
        imp2 = x(15);
        imp3 = x(16);
        imp4 = x(17);
        imp5 = x(18);
        imp6 = x(19);
        imp7 = x(20);
        imp8 = x(21);
        imp9 = x(22);
        
        ost_enrollment = x(23);
        nsp_enrollment = x(24);
        r_inc_up = x(25);
        start_year = x(26);
        
        r_F0F1 = x(27);
        r_F1F2 = x(28);
        r_F2F3 = x(29);
        r_F3F4  = x(30);
        r_F0F1_PWID  = x(31);
        r_F1F2_PWID  = x(32);
        r_F2F3_PWID  = x(33);
        r_F3F4_PWID = x(34);
        
    end



end
