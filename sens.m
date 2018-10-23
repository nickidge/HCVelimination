clear all

sens_base
sens_lowercov
sens_highercov
sens_highersensitivity
sens_annualtest
sens_90diagnosedtreated
sens_GPscreening

user=extractBetween(pwd,"Users\","\");
drive=extractBefore(pwd,":");
Table_base = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_base"),'paper2_text'));
Table_lowercov = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_lowercov"),'paper2_text'));
Table_highercov = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_highercov"),'paper2_text'));
Table_highersensitivity = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_highersensitivity"),'paper2_text'));
Table_annualtest = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_annualtest"),'paper2_text'));
Table_90diagnosedtreated = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_90diagnosedtreated"),'paper2_text'));
Table_GPscreening = struct2array(load(strcat(drive,":\Users\",user,"\Desktop\Matlab Sims\Tanzania\sens_GPscreening"),'paper2_text'));


sens_table = [Table_base([1,35:44],:);...
    Table_lowercov([1,35:44],:);...
    Table_highercov([1,35:44],:);...
    Table_highersensitivity([1,35:44],:);...
    Table_annualtest([1,35:44],:);...
    Table_90diagnosedtreated([1,35:44],:);...
    Table_GPscreening([1,35:44],:)];

    