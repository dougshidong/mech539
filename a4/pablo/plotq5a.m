%% a
close all;
clear;
let = 'a';
st = 4;
en = 26;
in = 2;
let = 'c';
st = 2;
en = 8;
in = 1;

nline = (en - st)/in + 1;

set(0,'DefaultTextInterpreter','Latex');

h = figure('Position', [0, 0, 750, 500]);

hold on;

linecolors = hsv(nline);
i=0;
for alfa = st:in:en
    i=i+1;
    fname = sprintf('./q5results/q5%cresult_%d.mat', let, alfa);
    load(fname);
    plot(z(1:nbp+1,1),Cpj,'-','color',linecolors(i,:));
    legend_str{i} = sprintf('X = %d', alfa);
end
grid on
sq = plot(NaN, NaN, 'sk');
i=i+1;
legend_str{i} = sprintf('Transition');
sq = plot(NaN, NaN, 'vk');
i=i+1;
legend_str{i} = sprintf('Laminar Separation');
sq = plot(NaN, NaN, '*k');
i=i+1;
legend_str{i} = sprintf('Turbulent Separation');

columnlegend(2,legend_str);

i=0;
for alfa = st:in:en
    i=i+1;
    fname = sprintf('./q5results/q5%cresult_%d.mat', let, alfa);
    load(fname);
    
    if resup(4) == 1
        [minvalue, minloc] = min(abs(resup(5)/100 - z(1:nbp/2)));
        plot(z(minloc),Cpj(minloc),'s','color',linecolors(i,:));
    end
    if reslo(4) == 1
        [minvalue, minloc] = min(abs(reslo(5)/100 - z(nbp/2+2:nbp+1)));
        minloc = minloc + nbp/2 + 1;
        plot(z(minloc),Cpj(minloc),'s','color',linecolors(i,:));
    end
    if resup(4) == 2
        [minvalue, minloc] = min(abs(resup(5)/100 - z(1:nbp/2)));
        plot(z(minloc),Cpj(minloc),'v','color',linecolors(i,:));
    end
    if reslo(4) == 2
        [minvalue, minloc] = min(abs(reslo(5)/100 - z(nbp/2+2:nbp+1)));
        minloc = minloc + nbp/2 + 1;
        plot(z(minloc),Cpj(minloc),'v','color',linecolors(i,:));
    end
    
    if resup(6) ~= 0
        [minvalue, minloc] = min(abs(resup(6)/100 - z(1:nbp/2)));
        plot(z(minloc),Cpj(minloc),'*','color',linecolors(i,:));
    end

    if reslo(6) ~= 0
        [minvalue, minloc] = min(abs(reslo(6)/100 - z(nbp/2+2:nbp+1)));
        minloc = minloc + nbp/2 + 1;
        plot(z(minloc),Cpj(minloc),'*','color',linecolors(i,:));
    end
end

set(gca,'YDir','Reverse');
title('Pressure Distribution for Varying Camber NACA4X12');
xlabel('x')
xlim([0,1])
ylabel('$C_p$')
expfname = sprintf('../report/figures/q5%cpressure.pdf', let);
export_fig (expfname)
hold off

