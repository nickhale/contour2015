clear all


x = .99;
x = .999;
x = .9999;
x = .99999;
x = .999999;

ell2 = pi^2;
nn = [1:100 200];

LW = 'LineWidth'; lw = 2; cols = get(0, 'DefaultAxesColorOrder');

%% Normal method

for n = nn
    [z, w, h] = nodesandweights(n, x, ell2);
    u = 0;
    for k = 1:n
        v = 1./z(k).*(1-1./cosh(.5*sqrt(z(k))));
        u = u + v*w(k);
    end
    U(n) = imag(u);
    hh(n) = h;
end

err = abs(U - U(end)).';
semilogy(nn, err(nn), '-', nn, exp(-pi^2./hh(nn)), '--', 'color', cols(1,:), LW, lw);
hold on

%% x --> 1 method

for n = nn
    [z, w] = nodesandweights2(n, x, ell2);
    u = 0;
    for k = 1:n
        v = 1./z(k).*(1 - 1./cosh(.5*sqrt(z(k))));
        u = u + (v - 1./(z(k)-pi^2))*w(k);
    end
    U(n) = imag(u);
end

err2 = abs(U - U(end)).';
semilogy(nn, err2(nn), '-',  nn, exp(-pi*sqrt(nn)), '--', 'color', cols(2,:), LW, lw);

%%
% Make pretty legends and print to file.
nn2 = 1:5:100;
MFC = 'MarkerFaceColor'; MS = 'MarkerSize';

ms = 6; lw = 2;
plot(nn2, err(nn2), 'o', MS, ms, 'color', cols(1,:), MFC, 'w', LW, 2)
h1 = plot(NaN, NaN, '-o', MS, ms, 'color', cols(1,:), MFC, 'w', LW, 2);
plot(nn2, err2(nn2), 'o', MS, ms, 'color', cols(2,:), MFC, cols(2,:), LW, 2)
h2 = plot(NaN, NaN, '-o', MS, ms, 'color', cols(2,:), MFC, cols(2,:), LW, 2);
hold off
if ( x == .99 )
    legend([h1 h2], 'h_\ast from (2.12)', 'h_\ast from (5.3)')
end
xlim([0 100]), ylim([1e-15, 1e0])
grid on
set(gca, 'fontsize', 16), shg
str = num2str(x);
if ( length(str) > 1 )
    str = str(3:end);
end
print('-depsc2', ['../paper/figures/x1lim' str])






