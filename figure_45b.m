NN = 1:20;
err = [];

% Solve the particle in a box problem.
for N = NN
    disp(N)
    b = .1;                                 % Aspect ratio of box.
    x = .5;                                 % Evaluate in centre of box.
    m = 20;                                 % Space discretization size.
    [~, L] = cheb2bc(m, [0 1 0 ; 1 0 0]);   % Neumann condition one side
                                            %   and Dirichlet at the other.
    A = (4/b)^2*L;                          % Scaling.
    U1 = 2*ones(m - 1);                     % RHS.
    [S, A] = schur(A, 'complex');           % Compute Schur factorisation only once.
    U1 = S'*U1*S;                           % Transform the RHS.
    ell2 = 2*(pi/b)^2;                      % Smallest eigenvalue of operator.
    [z, w, hh(N)] = contourIntegral(N, x, ell2);   % Quadrature nodes and weights.
    % Solve resolvent systems:
    U = 0; I = eye(m-1);
    for k = 1:N
        res = lyap( z(k)/2*I - A, z(k)/2*I - A', -U1 );
        U = U + imag(w(k)*res);
    end
    % Undo the Schur factorization:
    U = S*U*S';
    uApprox = U(1);                          % Computed solution.
    uExact = 7.2988176570485260889e-10;     % Know solution.
%     uExact = 1/3;
    relErr(N) = abs((uApprox - uExact)/uExact); % Relative Error.
end

z1 = z;

%%
ca
LW = 'LineWidth'; lw = 2; MFC = 'MarkerFaceColor';
semilogy(NN, relErr(NN), '-o', LW, 2, MFC, 'w'), shg, grid on
hold on
semilogy(NN, exp(-pi^2./hh(NN))/uExact, '--k', LW, 2);
ylim([1e-13, 1e5])
set(gca, 'FontSize', 16);
pause(1)
% print -depsc2 ~/Dropbox/work/2015/contour/paper/figures/piab_conv

return

%% COMPARE WITH HHT

NN = 5:1:80;

m = 20;                               % m = dimension of DM+2
[t, L] = cheb2bc(m, [0 1 0; 1 0 0]);  % Dirichlet-condition one side,
L = (4/b)^2*L;                        % Scaling.

f = @(z) sinh(x*sqrt(-z))./sinh(sqrt(-z)); % change this for another function f
I = eye(size(L));
U0 = 2*ones(size(L));

kronL = kron(L, I) + kron(I, L);
e = real(eig(kronL));
eR = max(e); eL = min(e);
% eR = 0; % Do not approx min eval.
ell2 = -eR + pi^2; L2 = -eL + pi^2;

k = (sqrt(L2/ell2)-1)/(sqrt(L2/ell2)+1);
[K, Kp] = ellipkkp(-log(k)/pi);
counter = 0;
for N = NN
    disp(N)
    counter = counter + 1;
    t = .5i*Kp - K + (.5:N)'*2*K/N;
    [u, cn, dn] = ellipjc(t, -log(k)/pi);
    z = -sqrt(ell2*L2)*((1/k+u)./(1/k-u))+pi^2;
    dzdt = -cn.*dn./(1/k-u).^2;
    w = -4*K*sqrt(ell2*L2)*f(z).*dzdt/(k*pi*N);
    
    idx = isnan(w);
    z(idx) = []; w(idx) = [];
    tmp = 0; tol = eps;
    
    U = 0;
    for j = length(w):-1:1
        dS = w(j)*lyap(L-z(j)/2*I, L'-z(j)/2*I, U0);
        if ( norm(imag(dS(1)), inf) > uExact*tol )
            tmp = tmp + 1;
            U = U + dS(:);
        else
            w(j) = [];
            z(j) = [];
        end
        
    end
    U = imag(U);
    uapp =  U(1);
    reler = abs((uapp-uExact)/uExact);
    err(counter) = reler;
    NN(counter) = tmp;
    
end
z2 = z;
semilogy(NN, err, '-', 'linewidth', 2), shg, grid on

%%

figure(2)
plot(z1, 'o'); hold on
plot(z2', '.r'); hold off
xlim([-2000, 1000])