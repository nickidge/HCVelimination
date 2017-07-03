function apply_cascade_rates(scenario, mult)
global c1_acutePWID c1_chronicPWID c2_PWID c3_PWID c4_PWID c5_PWID c6_PWID ...
    c1_acuteformerPWID c1_chronicformerPWID c2_formerPWID c3_formerPWID c4_formerPWID c5_formerPWID c6_formerPWID ...
    c1_acute c1_chronic c2 c3 c4 c5 c6 care ...
    c1_chronicPWID_base c2_PWID_base c3_PWID_base c4_PWID_base c5_PWID_base c1_chronicformerPWID_base c2_formerPWID_base c3_formerPWID_base c4_formerPWID_base c5_formerPWID_base c1_chronic_base c2_base c3_base c4_base c5_base ...
    filename

if strcmp(filename,'C:\Users\Nick\Desktop\Matlab Sims\Cascade\doubletime')==1
    specialist_fibroscan_rate = 0.25 * 1/(2*65/365);
else
    specialist_fibroscan_rate = 0.25 * 1/(65/365);
end
primary_fibroscan_rate = 0.5 * 1/(65/365);% 52; % appoximating instantly

if strcmp(scenario,'old') == 1 % Current ratio of GP / community
    
    c3_PWID = c3_PWID_base;
    c3_formerPWID = c3_formerPWID_base;
    c3 = c3_base;
    c4_PWID = c4_PWID_base;
    c4_formerPWID = c4_formerPWID_base;
    c4 = c4_base;

end

if strcmp(scenario,'initial') == 1 % Current ratio of GP / community
    
    c3_PWID0 = c3_PWID_base;
    c3_formerPWID0 = c3_formerPWID_base;
    c30 = c3_base;
    c4_PWID0 = c4_PWID_base;
    c4_formerPWID0 = c4_formerPWID_base;
    c40 = c4_base;
    
    c3_PWID1 = 100; % remove the need to wait for genotype
    c3_formerPWID1 = 100;
    c31 = 100;
    c4_PWID1 = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    c4_formerPWID1 = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    c41 = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    
    
    c3_PWID = c3_PWID0 + mult * (c3_PWID1 - c3_PWID0);
    c3_formerPWID = c3_formerPWID0 + mult * (c3_formerPWID1 - c3_formerPWID0);
    c3 = c30 + mult * (c31 - c30);
    c4_PWID = c4_PWID0 + mult * (c4_PWID1 - c4_PWID0);
    c4_formerPWID = c4_formerPWID0 + mult * (c4_formerPWID1 - c4_formerPWID0);
    c4 = c40 + mult * (c41 - c40);
    
end
if strcmp(scenario,'current') == 1 % Current ratio of GP / community
    
    %c4_PWID = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    %c4_formerPWID = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    %c4 = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    
end

if strcmp(scenario,'scale-up') == 1 || strcmp(scenario,'rapidRNA') == 1 % scales up from initial ratio of hospticat to community to final ratio of hospital to community


    c4_PWID0 = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    c4_formerPWID0 = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    c40 = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    
    c4_PWID1 = (1-care(2))* specialist_fibroscan_rate + care(2)* primary_fibroscan_rate;
    c4_formerPWID1 = (1-care(2))* specialist_fibroscan_rate + care(2)* primary_fibroscan_rate;
    c41 = (1-care(2))* specialist_fibroscan_rate + care(2)* primary_fibroscan_rate;
    
    c4_PWID = c4_PWID0 + mult * (c4_PWID1 - c4_PWID0);
    c4_formerPWID = c4_formerPWID0 + mult * (c4_formerPWID1 - c4_formerPWID0);
    c4 = c40 + mult * (c41 - c40);
end


if strcmp(scenario,'APRI') == 1 || strcmp(scenario,'WHO') == 1% genotype and liver disease done together


    c4_PWID0 = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    c4_formerPWID0 = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    c40 = (1-care(1))* specialist_fibroscan_rate + care(1)* primary_fibroscan_rate;
    
    c4_PWID1 = (1-care(2))* specialist_fibroscan_rate + care(2)* primary_fibroscan_rate;
    c4_formerPWID1 = (1-care(2))* specialist_fibroscan_rate + care(2)* primary_fibroscan_rate;
    c41 = (1-care(2))* specialist_fibroscan_rate + care(2)* primary_fibroscan_rate;
    
    c4_PWID = c4_PWID0 + mult * (c4_PWID1 - c4_PWID0);
    c4_formerPWID = c4_formerPWID0 + mult * (c4_formerPWID1 - c4_formerPWID0);
    c4 = c40 + mult * (c41 - c40);

end


end

% 
% if strcmp(scenario,'current') == 1 % Current ratio of GP / community
%     
%     c4_PWID = (1-care(1))* 1/(65/365) + care(1)* 1/(1/365);
%     c4_formerPWID = (1-care(1))* 1/(65/365) + care(1)* 1/(1/365);
%     c4 = (1-care(1))* 1/(65/365) + care(1)* 1/(1/365);
%     
% end
% 
% if strcmp(scenario,'fast') == 1 
%     c1_chronicPWID = 1/0.5;
%     c2_PWID = 10;
%     c3_PWID = 10;
%     c4_PWID = 10;
%     c5_PWID = 10;
%     c1_chronicformerPWID = 1/0.5;
%     c2_formerPWID = 10;
%     c3_formerPWID = 10;
%     c4_formerPWID = 10;
%     c5_formerPWID = 10;
%     c1_chronic = 1/0.5;
%     c2 = 10;
%     c3 = 10;
%     c4 = 10;
%     c5 = 10;
% end
% 
% if strcmp(scenario,'90-90-90') == 1
%     c1_chronicPWID = 1.3126;
%     c2_PWID = 2.0614 ;
%     c3_PWID = 2.0121;
%     c4_PWID = 1.5607;
%     c1_chronicformerPWID = 0.9363;
%     c2_formerPWID = 2.3652 ;
%     c3_formerPWID = 2.5831;
%     c4_formerPWID = 2.8806 ;
%     c1_chronic = 0.3335;
%     c2 = 0.1449;
%     c3 = 1.6845;
%     c4 = 2.3063;
% end
% 
% if strcmp(scenario,'nodeaths') == 1
%     c1_chronicPWID = 1/0.5;
%     c2_PWID = 10;
%     c3_PWID = 10;
%     c4_PWID = 10;
%     c1_chronicformerPWID = 1/0.5;
%     c2_formerPWID = 10;
%     c3_formerPWID = 10;
%     c4_formerPWID = 10;
%     c1_chronic = 1/0.5;
%     c2 = 10;
%     c3 = 10;
%     c4 = 10;
% end
