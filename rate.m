function ad = rate(r_XY_PWID,r_XY) % sub-function that distributes factors for PWID and everyone else
global num_pops num_cascade num_age num_intervention num_engagement num_region

ad = reshape(repmat([r_XY_PWID;r_XY;r_XY],num_cascade,num_age,num_intervention,num_engagement,num_region),num_pops,num_cascade,num_age,num_intervention, num_engagement,num_region);

end