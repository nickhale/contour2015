function test_bounds2

xx = .1:.2:.9;
ee = logspace(-6, -1, 21);
ell = pi;

    function f = F(t)
        f = abs( E(x, z(t)) .* R(z(t)) .* zp(t) );
        f(isnan(f)) = 0;
    end

% figure(1)
% for x = xx
%     z = @(t) .5*(pi^2-ell^2) + .5i*(pi^2+ell^2)*sinh(t);
%     zp = @(t) .5i*(pi^2+ell^2)*cosh(t);
%     Ip = [];
%     for k = 1:numel(ee)
%         k
%         Ip(k) = 2*integral(@(t) F(t+1i*(pi/2 - ee(k))), 0 , inf);
%     end
%     loglog(ee, Ip./log(1./ee), 'LineWidth', 2); hold on
% end
% 
% hold off, grid on
% set(gca, 'fontsize', 14)
% axis([1e-7, 1e0 1e-1, 1e1])

% print -depsc2 ../paper/figures/Mp

%%

figure(2)
for x = xx
    z = @(t) .5*(pi^2-ell^2) + .5i*(pi^2+ell^2)*sinh(t);
    zp = @(t) .5i*(pi^2+ell^2)*cosh(t);
    Im = [];
    for k = 1:numel(ee)
        k
        Im(k) = 2*integral(@(t) F(t-1i*(pi/2 - ee(k))), 0 , inf);
    end
    loglog(ee, Im./log(1./ee), '-', 'LineWidth', 2);  hold on
end

hold off, grid on
set(gca, 'fontsize', 14)
axis([1e-7, 1e0 1e0, 1e1])

legend('x = 0.1', 'x = 0.3', 'x = 0.5', 'x = 0.7', 'x = 0.9', 'Location', 'SW')

% print -depsc2 ../paper/figures/Mm

end

function r = R(z)
y = linspace(0, .5, 101);
v = @(y, z) 1./z .* (1 - cosh((.5-y)*sqrt(z))./cosh(.5*sqrt(z)) );
r = zeros(size(z));
for k = 1:numel(z)
    r(k) = norm(v(y, z(k)), inf);
end
end
