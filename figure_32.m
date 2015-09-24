% Error as a function of number of terms

clc, close all, clear all

x = .5;
x = .7;
% x = .95;

%%

NN = 1:52;
NN2 = 1:2:51;
ell2 = pi^2;
y = linspace(0, 1, 101);

%% 'Exact' solution (many terms)

sol = 0;
for k = 1:2:1000
    ds = 2/(pi*k)*(1-(-1)^k)*sinh(k*pi*x)*sin(k*pi*y)/sinh(k*pi);
    ds(isnan(ds)) = 0;
    sol = sol + ds;
    if ( ~any(ds(:)) ), break, end
end

%% Contour series:

l = 0;
for N = NN
    l = l + 1;
    [z, w, h] = nodesandweights(N, x, ell2);
    u = 0*y; 
    for k = 1:N
        zk = z(k);
        v = 1./zk.*(1 - cosh((.5-y)*sqrt(zk))./cosh(.5*sqrt(zk)));
        v(isinf(v) | isinf(v)) = 0;
        u = u + imag(v*w(k));
    end    
    err_con(l) = norm(u - sol, inf);
end

%% Classic series:

l = 0; s = 0;
for k = NN2
    l = l + 1;
    s = s + 2/(pi*k)*(1-(-1)^k)*sinh(k*pi*x)*sin(k*pi*y)/sinh(k*pi);
    err_ser(l) = norm(s - sol, inf);
end
    
%%

MFC = 'MarkerFaceColor';
LW = 'LineWidth'; lw = 3;
cols = get(0, 'DefaultAxesColorOrder');
blue = cols(1,:);
red  = cols(2,:);
MS = 'MarkerSize'; ms = 7;

semilogy(NN, err_con, '-', LW, lw, MFC, 'w', 'color', blue, MS, ms); hold on
h2 = semilogy(NN2, err_ser, '-o' ,LW, lw, MFC, red, 'color', red, MS, ms);

h1 = plot([NaN NaN], [NaN NaN],  '-o', 'color', blue, MFC, 'w', LW, lw, MS, ms);

semilogy(NN2+1, err_con(NN2+1), 'o', LW, lw, MFC, 'w', 'color', blue, MS, ms);

grid on
axis([0 40 1e-15, 1])
if ( x == .5 )
    legend([h1 h2], 'Contour integral', 'Classic series')
end
set(gca, 'fontsize', 16);

eval(['print -depsc2 ../paper/figures/contourVsClassic_x=' num2str(100*x) '.eps'])

