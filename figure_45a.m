% close all, clc, clear all
ca
LW = 'LineWidth'; lw = 3;

% mm = [3 3:31];
mm = unique(ceil(logspace(log10(3), log10(200), 50)));
% mm = 19;

b = .1;                                 % Aspect ratio of box.
x = .5;                                 % Evaluate in centre of box.
N = 20;                                 % N = number of resolvent solves.

trueSoln = 7.298817657294230e-10;
loopnum = @(m) max(ceil((200-m.^1.5)/20), 5);
loopnum = @(m) 1 + 0*m;

%%
t1 = zeros(mm(end), 1);
for m = [mm(1) mm]
    disp(m)
    tic
    for loop = 1:loopnum(m)
        [~, L] = cheb2bc(m, [0 1 0 ; 1 0 0]);   % Neumann condition one side
                                                %   and Dirichlet at the other.
        L = (4/b)^2*L;                          % Scaling.
        U1 = 2*ones(m - 1); U1 = U1(:);         % RHS.
        ell2 = (pi/b)^2;                        % Smallest eigenvalue of operator.
        [z, w] = contourIntegral(N, x, ell2);   % Quadrature nodes and weights.
        % Solve resolvent systems:
        U = 0; I = eye(m-1, m-1);
        A = kron(L, I) + kron(I, L);
        II = eye((m-1)^2, (m-1)^2);
        for k = 1:N
            res = (z(k)*II - A)\U1;
            U = U + imag(w(k)*res);
        end
        uApprox = U(1);                          % Computed solution.
    end
    tmp = toc;
    t1(m) = tmp/loopnum(m);
    if ( t1(m) > 1 ), break, end
end

%%
t2 = zeros(mm(end), 1);
for m = [mm(1) mm]
    disp(m)
    tic
    for loop = 1:loopnum(m)
        [~, L] = cheb2bc(m, [0 1 0 ; 1 0 0]);   % Neumann condition one side
                                                %   and Dirichlet at the other.
        L = (4/b)^2*L;                          % Scaling.
        U1 = 2*ones(m - 1); U1 = U1(:);         % RHS.
        ell2 = (pi/b)^2;                        % Smallest eigenvalue of operator.
        [z, w] = contourIntegral(N, x, ell2);   % Quadrature nodes and weights.
        % Solve resolvent systems:
        U = 0; I = eye(m-1, m-1);
        A = kron(L, I) + kron(I, L);
        [P, H] = hess(A);                       % Hessenberg factorization.
        U1 = P'*U1;
        II = eye((m-1)^2, (m-1)^2);
        for k = 1:N
            res = (z(k)*II - H)\U1;
            U = U + w(k)*res;
        end
        U = imag(P*U);
        uApprox = U(1);                          % Computed solution.
    end
    tmp = toc;
    t2(m) = tmp/loopnum(m);
    if ( t2(m) > 1 ), break, end
end


%%
% t2b = zeros(mm(end), 1);
% for m = mm
%     disp(m)
%     tic
%     for loop = 1:loopnum(m)
%         [~, L] = cheb2bc(m, [0 1 0 ; 1 0 0]);   % Neumann condition one side
%                                                 %   and Dirichlet at the other.
%         L = (4/b)^2*L;                          % Scaling.
%         U1 = 2*ones(m - 1); U1 = U1(:);         % RHS.
%         ell2 = (pi/b)^2;                        % Smallest eigenvalue of operator.
%         [z, w] = contourIntegral(N, x, ell2);   % Quadrature nodes and weights.
%         % Solve resolvent systems:
%         U = 0; I = eye(m-1, m-1);
%         A = kron(L, I) + kron(I, L);
%         [P, H] = schur(A, 'complex');           % Schur factorization.
%         U1 = P'*U1;
%         II = eye((m-1)^2, (m-1)^2);
%         for k = 1:N
%             res = (z(k)*II - H)\U1;
%             U = U + w(k)*res;
%         end
%         U = imag(P*U);
%         uApprox = U(1);                          % Computed solution.
%     end
%     tmp = toc;
%     t2b(m) = tmp/loopnum(m);
%     if ( t2b(m) > 1 ), break, end
% end


%%

t3 = zeros(mm(end), 1);
for m = [mm(1) mm]
    disp(m)
    tic
    for loop = 1:loopnum(m)
        [~, L] = cheb2bc(m, [0 1 0 ; 1 0 0]);   % Neumann condition one side
                                                %   and Dirichlet at the other.
        A = (4/b)^2*L;                          % Scaling.
        U1 = 2*ones(m - 1);                     % RHS.
        ell2 = (pi/b)^2;                        % Smallest eigenvalue of operator.
        [z, w] = contourIntegral(N, x, ell2);   % Quadrature nodes and weights.
        % Solve resolvent systems:
        U = 0; I = eye(m-1);
        for k = 1:N
            res = lyap( z(k)/2*I - A, z(k)/2*I - A', -U1 );
            U = U + imag(w(k)*res);
        end
        uApprox = U(1);                          % Computed solution.
    end
    tmp = toc;
    t3(m) = tmp/loopnum(m);
    if ( t3(m) > 1 ), break, end
end


%%

t4 = zeros(mm(end), 1);
for m = [mm(1) mm]
    disp(m)
    tic
    for loop = 1:loopnum(m)
        [~, L] = cheb2bc(m, [0 1 0 ; 1 0 0]);   % Neumann condition one side
                                                %   and Dirichlet at the other.
        A = (4/b)^2*L;                          % Scaling.
        U1 = 2*ones(m - 1);                     % RHS.
        I = eye(m-1);
        [S, A] = schur(A, 'complex');           % Compute Schur factorisation only once.
        U1 = S'*U1*S;                           % Transform the RHS.
        ell2 = (pi/b)^2;                        % Smallest eigenvalue of operator.
        [z, w] = contourIntegral(N, x, ell2);   % Quadrature nodes and weights.
        % Solve resolvent systems:
        U = 0; I = eye(m-1);
        for k = 1:N
            res = lyap( z(k)/2*I - A, z(k)/2*I - A', -U1 );
            U = U + imag(w(k)*res);
        end
        % Undo the Schur factorization:
        U = S*U*S';
        uApprox = U(1);                          % Computed solution.
    end
    tmp = toc;
    t4(m) = tmp/loopnum(m);
    if ( t4(m) > 1 ), break, end
end


%%
ca
cols = get(0, 'DefaultAxesColorOrder');
MFC = 'MarkerFaceColor'; MS = 'MarkerSize';

h1 = loglog(mm, t1(mm), LW, lw, 'color', cols(1,:)); shg, hold on
h2 = loglog(mm, t2(mm), ':', LW, lw, 'color', cols(1,:));
% h2b = loglog(mm, t2b(mm), '--', LW, lw); 
h3 = loglog(mm, t3(mm), LW, lw, 'color', cols(2,:)); 
h4 = loglog(mm, t4(mm), ':', LW, lw, 'color', cols(2,:));

% mm2 = ceil(logspace(log10(mm(1)), log10(mm(end)), 5));
mm2 = mm(1:5:end)
mm2(2) = 6; mm2(3) = 11;

h1b = loglog(NaN, NaN, '-s', MFC, 'w', MS, 8, 'color', cols(1,:), LW, 2);
loglog(mm2, t1(mm2), 's', MFC, 'w', MS, 8, 'color', cols(1,:), LW, 2);
h2b = loglog(NaN, NaN, '--s', MFC, 'w', MS, 8, 'color', cols(1,:), LW, 2);
loglog(mm2, t2(mm2), 's', MFC, 'w', MS, 8, 'color', cols(1,:), LW, 2);
h3b = loglog(NaN, NaN, '-^', MFC, 'w', MS, 8, 'color', cols(2,:), LW, 2);
loglog(mm2, t3(mm2), '^', MFC, 'w', MS, 8, 'color', cols(2,:), LW, 2);
h4b = loglog(NaN, NaN, '--^', MFC, 'w', MS, 8, 'color', cols(2,:), LW, 2);
loglog(mm2, t4(mm2), '^', MFC, 'w', MS, 8, 'color', cols(2,:), LW, 2);

legend([h1b h2b h3b h4b], '2D Kron', '2D Kron + Hessenberg', 'Lyap', 'Lyap + Schur', 'location', 'NW')
grid on
ylim([5e-4, 2e0])
xlim([1, 3e2])
set(gca, 'FontSize', 16);

%%
% pause(1)
% print -depsc2 ~/Dropbox/work/2015/contour/paper/figures/piab_time
% savefig ~/Dropbox/work/2015/contour/paper/figures/piab_time.fig

