function rho = SDCPupdate(SDCP_old,pseudo_new1,pseudo_new2,carrier_new1,carrier_old1,carrier_new2,carrier_old2,num_steps)

M = num_steps;

rho = 1/M*(pseudo_new2 - pseudo_new1) + (M-1)/M*(SDCP_old + carrier_new2 - carrier_old2 + carrier_old1 - carrier_new1);


end