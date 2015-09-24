LW = 'LineWidth';

s = .005;
grey =.9*[1 1 1];
fill([s 1-s 1-s s s], [s s 1-s 1-s s], grey ); hold on
plot([0 1 1 0 0], [0 0 1 1 0], 'k',  LW, 2);

s = .07;
plot([-2*s, 2*s], [-s, -s], 'k')
plot([-s, -s], [-2*s, 2*s], 'k')
plot(2*[s-.01 s s-.01], -[s+.01 s s-.01], 'k')
plot(-[s+.01 s s-.01], 2*[s-.01 s s-.01], 'k')

hold off

axis off equal square
axis([-.1 1.1, -.1 1.1])

% print -depsc2 domain1.eps