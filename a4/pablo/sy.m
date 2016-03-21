function cd = sy(thup,Hup,ueteup,thlo,Hlo,uetelo);

cd = 2*thup*(ueteup).^((Hup+5)/2) + 2*thlo*(uetelo).^((Hlo+5)/2);

cd = floor(round(10000*cd))/10000;
