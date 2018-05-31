global  delta alpha p_complete... 
    r_AF0 r_F0F1 r_F1F2 r_F2F3 r_F3F4 r_F0F1_PWID r_F1F2_PWID r_F2F3_PWID r_F3F4_PWID ...
    r_F4DC r_DCHCC r_F4HCC r_DCLT r_DCdeath r_HCCLT r_HCCdeath r_LTdeath1 r_LTdeath2 r_S4death r_LT1LT2...
    Q_sens q_svr q_svr_PWID q_treat 

load('calibration_test3','delta', 'alpha', 'p_complete',... 
    'r_AF0', 'r_F0F1', 'r_F1F2', 'r_F2F3', 'r_F3F4', 'r_F0F1_PWID', 'r_F1F2_PWID', 'r_F2F3_PWID', 'r_F3F4_PWID', ...
    'r_F4DC', 'r_DCHCC', 'r_F4HCC', 'r_DCLT', 'r_DCdeath', 'r_HCCLT', 'r_HCCdeath', 'r_LTdeath1', 'r_LTdeath2', 'r_S4death', 'r_LT1LT2',...
    'Q_sens', 'q_svr', 'q_svr_PWID', 'q_treat')
variance = 0.0001;

% alpha=random('Uniform',0.9,0.9);
% p_complete=random('Uniform',0.85,0.95);
% r_AF0 = 52/random(makedist('Triangular','a',r_AF0-variance*r_AF0,'b',r_AF0,'c',r_AF0+variance*r_AF0));
% r_F0F1 = random(makedist('Triangular','a',r_F0F1*(1-variance),'b',r_F0F1,'c',r_F0F1*(1+variance)));
r_F0F1 = -log(1-random(makedist('Triangular','a',1-exp(-r_F0F1*(1-variance)),'b',1-exp(-r_F0F1),'c',1-exp(-r_F0F1*(1+variance)))));
r_F1F2 = -log(1-random(makedist('Triangular','a',1-exp(-r_F1F2*(1-variance)),'b',1-exp(-r_F1F2),'c',1-exp(-r_F1F2*(1+variance)))));
r_F2F3 = -log(1-random(makedist('Triangular','a',1-exp(-r_F2F3*(1-variance)),'b',1-exp(-r_F2F3),'c',1-exp(-r_F2F3*(1+variance)))));
r_F3F4 = -log(1-random(makedist('Triangular','a',1-exp(-r_F3F4*(1-variance)),'b',1-exp(-r_F3F4),'c',1-exp(-r_F3F4*(1+variance)))));
r_F0F1_PWID = -log(1-random(makedist('Triangular','a',1-exp(-r_F0F1_PWID*(1-variance)),'b',1-exp(-r_F0F1_PWID),'c',1-exp(-r_F0F1_PWID*(1+variance)))));
r_F1F2_PWID = -log(1-random(makedist('Triangular','a',1-exp(-r_F1F2_PWID*(1-variance)),'b',1-exp(-r_F1F2_PWID),'c',1-exp(-r_F1F2_PWID*(1+variance)))));
r_F2F3_PWID = -log(1-random(makedist('Triangular','a',1-exp(-r_F2F3_PWID*(1-variance)),'b',1-exp(-r_F2F3_PWID),'c',1-exp(-r_F2F3_PWID*(1+variance)))));
r_F3F4_PWID = -log(1-random(makedist('Triangular','a',1-exp(-r_F3F4_PWID*(1-variance)),'b',1-exp(-r_F3F4_PWID),'c',1-exp(-r_F3F4_PWID*(1+variance)))));
r_F4DC = -log(1-random(makedist('Triangular','a',1-exp(-r_F4DC*(1-variance)),'b',1-exp(-r_F4DC),'c',1-exp(-r_F4DC*(1+variance)))));
r_DCHCC = -log(1-random(makedist('Triangular','a',1-exp(-r_DCHCC*(1-variance)),'b',1-exp(-r_DCHCC),'c',1-exp(-r_DCHCC*(1+variance)))));
r_F4HCC = -log(1-random(makedist('Triangular','a',1-exp(-r_F4HCC*(1-variance)),'b',1-exp(-r_F4HCC),'c',1-exp(-r_F4HCC*(1+variance)))));
r_HCCLT = -log(1-random(makedist('Triangular','a',1-exp(-r_HCCLT*(1-variance)),'b',1-exp(-r_HCCLT),'c',1-exp(-r_HCCLT*(1+variance)))));
r_DCLT = -log(1-random(makedist('Triangular','a',1-exp(-r_DCLT*(1-variance)),'b',1-exp(-r_DCLT),'c',1-exp(-r_DCLT*(1+variance)))));
r_DCdeath = -log(1-random(makedist('Triangular','a',1-exp(-r_DCdeath*(1-variance)),'b',1-exp(-r_DCdeath),'c',1-exp(-r_DCdeath*(1+variance)))));
r_HCCdeath = -log(1-random(makedist('Triangular','a',1-exp(-r_HCCdeath*(1-variance)),'b',1-exp(-r_HCCdeath),'c',1-exp(-r_HCCdeath*(1+variance)))));
r_LTdeath1 = -log(1-random(makedist('Triangular','a',1-exp(-r_LTdeath1 *(1-variance)),'b',1-exp(-r_LTdeath1 ),'c',1-exp(-r_LTdeath1 *(1+variance)))));
r_LTdeath2 = -log(1-random(makedist('Triangular','a',1-exp(-r_LTdeath2 *(1-variance)),'b',1-exp(-r_LTdeath2),'c',1-exp(-r_LTdeath2 *(1+variance)))));
r_S4death = -log(1-random(makedist('Triangular','a',1-exp(-r_S4death*(1-variance)),'b',1-exp(-r_S4death),'c',1-exp(-r_S4death*(1+variance)))));
r_LT1LT2=1;


% r_F0F1_PWID = random(makedist('Triangular','a',r_F0F1_PWID*(1-variance),'b',r_F0F1_PWID,'c',r_F0F1_PWID*(1+variance)));
% r_F1F2 = random(makedist('Triangular','a',r_F1F2*(1-variance),'b',r_F1F2,'c',r_F1F2*(1+variance)));
% r_F1F2_PWID = random(makedist('Triangular','a',r_F1F2_PWID*(1-variance),'b',r_F1F2_PWID,'c',r_F1F2_PWID*(1+variance)));
% r_F2F3 = random(makedist('Triangular','a',r_F2F3*(1-variance),'b',r_F2F3,'c',r_F2F3*(1+variance)));
% r_F2F3_PWID = random(makedist('Triangular','a',r_F2F3_PWID*(1-variance),'b',r_F2F3_PWID,'c',r_F2F3_PWID*(1+variance)));
% r_F3F4 = random(makedist('Triangular','a',r_F3F4*(1-variance),'b',r_F3F4,'c',r_F3F4*(1+variance)));
% r_F3F4_PWID = random(makedist('Triangular','a',r_F3F4_PWID*(1-variance),'b',r_F3F4_PWID,'c',r_F3F4_PWID*(1+variance)));
% r_F4DC = random(makedist('Triangular','a',r_F4DC*(1-variance),'b',r_F4DC,'c',r_F4DC*(1+variance)));
% r_DCHCC = random(makedist('Triangular','a',r_DCHCC*(1-variance),'b',r_DCHCC,'c',r_DCHCC*(1+variance)));
% r_F4HCC = random(makedist('Triangular','a',r_F4HCC*(1-variance),'b',r_F4HCC,'c',r_F4HCC*(1+variance)));
% r_HCCLT = random(makedist('Triangular','a',r_HCCLT*(1-variance),'b',r_HCCLT,'c',r_HCCLT*(1+variance)));
% r_DCLT = random(makedist('Triangular','a',r_DCLT*(1-variance),'b',r_DCLT,'c',r_DCLT*(1+variance)));
% r_DCdeath = random(makedist('Triangular','a',r_DCdeath*(1-variance),'b',r_DCdeath,'c',r_DCdeath*(1+variance)));
% r_HCCdeath = random(makedist('Triangular','a',r_HCCdeath*(1-variance),'b',r_HCCdeath,'c',r_HCCdeath*(1+variance)));
% r_LTdeath1 = random(makedist('Triangular','a',r_LTdeath1*(1-variance),'b',r_LTdeath1,'c',r_LTdeath1*(1+variance)));
% r_LTdeath2 = random(makedist('Triangular','a',r_LTdeath2*(1-variance),'b',r_LTdeath2,'c',r_LTdeath2*(1+variance)));
% r_S4death=-log(1-random(truncate(makedist('Normal',0.02,.001),0.01,0.03)));
% r_LT1LT2=1;
% % 
% r_F0F1=-log(1-random(truncate(makedist('Normal',0.106,.028),0.094,0.205)));
% r_F0F1_PWID=-log(1-random(truncate(makedist('Normal',0.116,.042),0.059,0.228)));
% r_F1F2=-log(1-random(truncate(makedist('Normal',0.074,.028),0.064,0.175)));
% r_F1F2_PWID=-log(1-random(truncate(makedist('Normal',0.085,.011),0.065,0.110)));
% r_F2F3=-log(1-random(truncate(makedist('Normal',0.106,.033),0.092,0.225)));
% r_F2F3_PWID=-log(1-random(truncate(makedist('Normal',0.085,.025),0.049,0.147)));
% r_F3F4=-log(1-random(truncate(makedist('Normal',0.105,.024),0.092,0.187)));
% r_F3F4_PWID=-log(1-random(truncate(makedist('Normal',0.13,.067),0.053,0.319)));
% r_F4DC=-log(1-random(truncate(makedist('Normal',0.037,.016),0.030,0.092)));
% r_DCHCC=-log(1-random(truncate(makedist('Normal',0.068,.015),0.041,0.099)));
% r_F4HCC=-log(1-random(truncate(makedist('Normal',0.01,.007),0.009,0.038)));
% r_HCCLT=-log(1-random(truncate(makedist('Normal',0.1,.033),0.050,0.180)));
% r_DCLT=-log(1-random(truncate(makedist('Normal',0.033,.008),0.017,0.049)));
% r_DCdeath=-log(1-random(truncate(makedist('Normal',0.138,.032),0.074,0.202)));
% r_HCCdeath=-log(1-random(truncate(makedist('Normal',0.605,.033),0.545,0.676)));
% r_LTdeath1=-log(1-random(truncate(makedist('Normal',0.169,.021),0.127,0.210)));
% r_LTdeath2=-log(1-random(truncate(makedist('Normal',0.034,.005),0.024,0.043)));
% r_S4death=-log(1-random(truncate(makedist('Normal',0.02,.001),0.01,0.03)));
% r_LT1LT2=1;
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

