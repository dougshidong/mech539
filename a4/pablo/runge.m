function ynp1 = runge(dx,y,Re,uei,duedsi,ueip1,duedsip1);

tsep = 0;

% 1st stage

yt(1) = y(1);
yt(2) = y(2);

H1 = yt(2);
H = HofH1(H1);

if H == 3 
   tsep = 1;
end;

if tsep ==0

  rtheta = Re*uei*yt(1);

  yp(1) = -(H+2)*yt(1)*duedsi./uei + 0.5*cfturb(rtheta,H);
  yp(2) = -H1*(duedsi./uei + yp(1)./yt(1))+0.0306*(H1-3).^(-0.6169)./yt(1);

  yt(1) = y(1) + dx*yp(1);
  yt(2) = y(2) + dx*yp(2);

  ys(1) = y(1) + 0.5*dx*yp(1);
  ys(2) = y(2) + 0.5*dx*yp(2);

  % 2nd stage

  H1 = yt(2);
  H = HofH1(H1);

  if H == 3 
    tsep = 1;
  end;

  if tsep ==0

     rtheta = Re*ueip1*yt(1);

     yp(1,1) = -(H+2)*yt(1)*duedsip1./ueip1 + 0.5*cfturb(rtheta,H);
     yp(2,1) = -H1*(duedsip1./ueip1 + yp(1)./yt(1))+0.0306*(H1-3).^(-0.6169)./yt(1);

     ynp1(1) = ys(1) + 0.5*dx*yp(1);
     ynp1(2) = ys(2) + 0.5*dx*yp(2);

  end;

end;

if tsep ==1
  ynp1(1) = -2; % so that H == 3 in the main
  ynp1(2) = 0;
end;
