function H = fH(lambda);

if lambda < 0

  if lambda==-0.14
    lambda=-0.139;
    disp('H(lambda) : Lambda = -0.14 -> Lambda = -0.139'); 
  end;

  H = 2.088 + 0.0731./(lambda+0.14);

elseif lambda >= 0 

  H = 2.61 - 3.75*lambda + 5.24*lambda.^2;

end;

