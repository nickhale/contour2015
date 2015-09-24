% Generate figures xstar=...
%
% Fig. 4.1. Relative error for PIAB problem when using x ∗ to determine h ∗ .
% Solid = numerical error (using spectral collocation with m = 32, 38, and 83,
% respecively). Dashed = theoretical error (T + D + + D − ). The error is the
% maximum error on the collocation nodes.

clc
close all
clear all

ell2 = pi^2;
lw = 3;

xstar = .5; m = 32;
xstar = .7; m = 38;
xstar = .95; m = 83;

f = ones(m-2, 1);
[y, A] = cheb2bc(m, [1 0 0 ; 1 0 0]); % } <-- DMSUITE
I = eye(m-2);
A = 4*A; y = (y+1)/2;
xx = linspace(0, 1, 101); xx([1, end]) = [];

%%

for N = 10:10:30
    
    l = 0; j = 0;
    for x = xx
        
        
        %% 'Exact' solution (many terms)
        sol = 0;
        for k = 1:2:1000
            ds = 2/(pi*k)*(1-(-1)^k)*sinh(k*pi*x)*sin(k*pi*y)/sinh(k*pi);
            ds(isnan(ds)) = 0;
            if ( ~any(ds(:)) )
                break
            end
            sol = sol + ds;
        end
        sol(end) = 0;
        
        j = j+1;
        l = l + 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optimal step size
        h = 2/N*wapr(sqrt(2)*pi^2*N/(1-xstar)/(sqrt(pi^2+ell2)));
        % Nodes and weights:
        mu = (ell2+pi^2)/2;
        theta = ((0:N-1)+0.5)*h;
        z = mu*(1+1i*sinh(theta)) - ell2;
        dz = 1i*mu*cosh(theta);
        w = h/pi*dz.*sin(x*sqrt(z))./sin(sqrt(z));
        w(isnan(w)) = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        u =  0*y;
        for k = 1:N
            zk = z(k);
            v = (zk*I - A)\f;
            u = u + imag(v*w(k));
        end
        
        tmp = (u - sol)./sol; tmp([1,end]) = 0;
        err(j) = norm(tmp, inf);
        
        errt(j) = exp(-(1-x)*sqrt(pi^2+ell2)*exp(N*h/2)/(2*sqrt(2)));
        errd(j) = exp(-pi^2/h);
        
    end
    
    cols = get(0, 'DefaultAxesColorOrder');
    h1 = semilogy(xx, err, '-', 'LineWidth', 3); hold on
    h2 = semilogy(xx, errt + 2*errd, '--', 'LineWidth', 3);
    if ( N == 30 )
        N = 50 ;
    end
    set(h2, 'color', cols(N/10,:));
    set(h1, 'color', cols(N/10,:));
    
end

%%

grid on
xlim([0 1])
ylim([1e-15, 1])
set(gca, 'fontsize', 16);

pause(1)
eval(['print -depsc2 ../paper/figures/xstar=' num2str(100*xstar) '.eps'])

        