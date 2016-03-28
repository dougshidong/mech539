close all;

set(0,'DefaultTextInterpreter','Latex');

h = figure('Position', [0, 0, 750, 500]);
hold on;
for alfa = 0:10:10
    fname = sprintf('./q2results/q2result_inviscid%d.mat', alfa);
    load(fname);
    plot(z(1:nbp+1,1),Cpj,'-');
    
    fname = sprintf('./q2results/q2result_viscous%d.mat', alfa);
    load(fname);
    plot(z(1:nbp+1,1),Cpj,'-');
    
    fname = sprintf('./q2results/q2exp%d.mat', alfa);
    load(fname);
    plot(z(1:nbp+1,1),Cpj,'s');
end
grid on
legend('Inviscid \alpha = 0','Viscous \alpha = 0', ...
       'Experimental \alpha = 0',...
       'Inviscid \alpha = 10','Viscous \alpha = 10',...
       'Experimental \alpha = 10');
set(gca,'YDir','Reverse');
title('Pressure Validation');
xlabel('x')
ylabel('$C_p$')
export_fig ('../report/figures/q2pressure.pdf')
hold off