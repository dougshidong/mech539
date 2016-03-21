function z  = ellipse(ratio,nbpo2);

beta = 0:pi/nbpo2:2*pi;
el =0.5*(1+ real(exp(i*beta))+i*ratio*imag(exp(i*beta))).';

z = [real(el) imag(el)];
