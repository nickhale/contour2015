LW   = 'LineWidth';
MFC  = 'MarkerFaceColor';
cols = get(0, 'DefaultAxesColorOrder');
blue = cols(1,:);
red  = cols(2,:);
grey = .7*[1 1 1];

j = 1:5;
t = linspace(0,20, 1000);

plot([-90, 1e6], [0 0], '--', 'color', grey); hold on
plot([-100 -90], [0 0 ], ':', 'color', grey)
plot([0, 0], [-1e6, 1e6], '--', 'color', grey);

plot([-95 -pi^2], [0 0 ], '-', 'color', blue, LW, 2)
plot([-100 -95], [0 0 ], ':', 'color', blue, LW, 2)

plot(-j.^2*pi^2, zeros(size(j)), 'ok', MFC, blue, ...
    'color', blue, 'markersize', 8, LW, 1.5);
plot([pi^2, 45], [0 0 ], '-', 'color', red, LW, 2)
plot([45 50], [0 0 ], ':', 'color', red, LW, 2)
plot(j.^2*pi.^2, zeros(size(j)), 'ok', MFC, 'w', ...
    'color', red, 'markersize', 8, LW, 1.5);

plot(-t.^2, t, '--k', LW, 3);
plot(-t.^2, -t, '--k', LW, 3);

t1 = 4.125;
x1 = -t1.^2;
y1 = t1;
plot(x1+[0 0 3], y1+[-.4 0 0], '-k', LW, 3);

t1 = -3.5;
x1 = -t1.^2;
y1 = t1;
plot(x1+[0 0 -3], y1+[-.4 0 0], '-k', LW, 3);

hold off

axis off equal square
axis([-100, 50 -10 10])
shg

pause(1)
print -depsc2 ../paper/figures/fig1.eps