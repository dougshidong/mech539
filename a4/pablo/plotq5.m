%% a
close all;
clear;

set(0,'DefaultTextInterpreter','Latex');

h = figure('Position', [0, 0, 750, 500]);

hold on;

linecolors = hsv(13);
i=0;
for alfa = 4:2:26
    i=i+1;
    fname = sprintf('./q5results/q5aresult_%d.mat', alfa);
    load(fname);
    plot(z(1:nbp+1,1),Cpj,'-','color',linecolors(i,:));
    legend_str{i} = sprintf('XX = %d', alfa);
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

columnlegend(2,legend_str,'SouthEast');

i=0;
for alfa = 4:2:26
    i=i+1;
    fname = sprintf('./q5results/q5aresult_%d.mat', alfa);
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
title('Pressure Distribution for Varying Thickness NACA00XX');
xlabel('x')
ylabel('$C_p$')
export_fig ('../report/figures/q5bpressure.pdf')
hold off