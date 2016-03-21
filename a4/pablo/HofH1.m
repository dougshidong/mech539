function H=HofH1(H1);

if H1 <= 3.32
  H = 3;
elseif H1 < 5.3
  H = 0.6778 + 1.1536*(H1-3.3).^(-0.326);
else
  H = 1.1 + 0.86*(H1-3.3).^(-0.777);
end

