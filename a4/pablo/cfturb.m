function cf = cfturb(rtheta,H);

cf = 0.246*(10.^(-0.678*H))*rtheta.^(-0.268);

