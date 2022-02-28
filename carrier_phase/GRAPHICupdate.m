function rho = GRAPHICupdate(GRAPHIC_old,pseudo_new,pseudo_old,carrier_new,carrier_old,num_steps)

M = num_steps;

rho = (M+1)/(2*M)*pseudo_new - (M-1)/(2*M)*pseudo_old + (M-1)/M*(GRAPHIC_old + 0.5*(carrier_new-carrier_old));


end