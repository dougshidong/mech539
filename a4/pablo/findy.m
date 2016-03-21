function y = findy(z,x,side);

nbpo2 = (size(z,1)-1)/2;

if side ==1  % upper side

  inda = find(z(1:nbpo2+1,1)<x);
  ind = min(inda);
  y = z(ind,2) + (x-z(ind,1))/(z(ind-1,1)-z(ind,1))*(z(ind-1,2)-z(ind,2));

else          % lower side

  inda = find(z(nbpo2+1:2*nbpo2+1,1)<x);
  ind = max(inda) + nbpo2+1;
  y = z(ind,2) + (x-z(ind,1))/(z(ind+1,1)-z(ind,1))*(z(ind+1,2)-z(ind,2));

end;
