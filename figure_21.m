clc, close all

d = pi/6;
ell = pi;
a = 0;

LW  = 'LineWidth'; lw = 3;
MFC = 'MarkerFaceColor';
cols = get(0, 'DefaultAxesColorOrder');
blue = cols(1,:);
red  = cols(2,:);
grey =.9*[1 1 1];

w = linspace(-10, 10, 1000);
wp = w + 1i*pi/2;
wm = w - 1i*pi/2;

k = 0;
N = 7;
for s = linspace(-pi/2, pi/2, N);
    k = k + 1;
    ww{k} = w + 1i*s + 1i*eps;
end
ww([1, (N+1)/2, end]) = [];
N = N - 3;

%%

figure(1)
fill([-10 10 10 -10 -10], pi/2*[-1 -1 1 1 -1], grey, ...
    'edgecolor', 'none'); hold on
for k = 1:N
    plot(ww{k}, '--k', LW, 1); hold on
end
plot(wp, LW, lw, 'color', blue); hold on
plot(wm, LW, lw, 'color', red); 
plot(w + 1i*eps, '--k', LW, lw)
plot(0, pi/2, 'o', 'color', blue, MFC, blue);
plot(0, -pi/2, 'o', 'color', red, MFC, 'w');
trim = 5;
plot(0, -trim, '.w', 'markerfacecolor', 'none', 'markersize', eps);
plot(0,  trim, '.w', 'markerfacecolor', 'none', 'markersize', eps);
axis equal off
axis(trim*[-1 1 -1 1])
set(gca, 'xtick', [],  'ytick', [])
hold off

print('-depsc2', ['../paper/figures/conformalmap_strip_a=', num2str(10*a)])

%%

figure(2)
z = @(w) .5*(ell^2-pi^2) + .5*1i*(pi^2+ell^2)*sinh(w);
c = 20;
b = (20-ell^2)*tan(d);
fill([-c c c -c -c ], [-c -c c c -c], grey, ...
    'edgecolor', 'none'); hold on
for k = 1:N
    plot(z(ww{k}), '--k', LW, 1); hold on
end
plot(z(wp), LW, lw, 'color', blue); hold on
plot(z(wm), LW, lw, 'color', red); 
plot(z(w), '--k', LW, lw)
plot(-ell^2, 0, 'o', 'color', blue, MFC, c1);
plot(pi^2,   0, 'o', 'color', red, MFC, 'w');
axis equal off
axis([-20 20 -20 20])
set(gca, 'xtick', [], 'ytick', [], 'box', 'off')
hold off

print('-depsc2', ['../paper/figures/conformalmap_elliptic_a=', num2str(10*a)])

%%

alignfigs

