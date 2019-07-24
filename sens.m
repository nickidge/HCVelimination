clear all

sens_base
sens_lowercov
sens_highercov
sens_highersensitivity
sens_annualtest
sens_UHC0
sens_UHC1
sens_50diagnosedtreated
sens_90diagnosedtreated
sens_GPscreening

user=extractBetween(pwd,"Users\","\");
drive=extractBefore(pwd,":");
Table_base = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_base"),'paper2_text'));
Table_lowercov = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_lowercov"),'paper2_text'));
Table_highercov = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_highercov"),'paper2_text'));
Table_highersensitivity = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_highersensitivity"),'paper2_text'));
Table_annualtest = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_annualtest"),'paper2_text'));
Table_UHC0 = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_UHC0"),'paper2_text'));
Table_UHC1 = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_UHC0"),'paper2_text'));
Table_50diagnosedtreated = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_50diagnosedtreated"),'paper2_text'));
Table_90diagnosedtreated = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_90diagnosedtreated"),'paper2_text'));
Table_GPscreening = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_GPscreening"),'paper2_text'));


sens_table = [Table_base([1,41:52],:);...
    Table_lowercov([1,41:52],:);...
    Table_highercov([1,41:52],:);...
    Table_highersensitivity([1,41:52],:);...
    Table_annualtest([1,41:52],:);...
    Table_UHC0([1,41:52],:);...
    Table_UHC1([1,41:52],:);...
    Table_50diagnosedtreated([1,41:52],:);...
    Table_90diagnosedtreated([1,41:52],:);...
    Table_GPscreening([1,41:52],:)];

    