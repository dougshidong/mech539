function L = fL(lambda);

if lambda < 0

  if lambda==-0.107
    lambda=-0.106;
    disp('l(lambda) : Lambda = -0.107 -> Lambda = -0.106');
  end;

  L = 0.22 + 1.402*lambda + (0.018*lambda)./(lambda+0.107);

elseif lambda >= 0 

  L = 0.22 + 1.57*lambda - 1.8*lambda.^2;

end;

