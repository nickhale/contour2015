function u = uexact(z)
% 'exact' solutions to particle in a box problem.

if ( z == .35 )
    u = 9.310504159555458e-13;       % z = .35;
elseif ( z == .5 )
    u = .72988176570485260889e-9;    % z = .5;
elseif ( z == .75)
    u = 4.864643935634531e-05;
elseif ( z == .9 )
    u = 0.038031567284990;
elseif ( z == .95 )
    u = 0.337420432289655;           % z = .95
elseif ( z == .97 )
    u = 0.7600263254;
elseif ( z == .99 )
    u = 1.5213585869;
elseif ( z == .991 )
    u = 1.567519810368169;
elseif ( z == .999 )
    u = 1.951131218785719;
    u = 1.951126560685;
    u = 1.951127946;
elseif ( z == .9999 )
    u = 1.995102638821978;
elseif ( z == .99999 )
    u = 1.999509233423;
elseif ( z == .999999 )
    u = 1.9999509094;
elseif ( z == .9999999 )
    u = 1.99999509073;
elseif ( z == 1 )
    u =  2;
else
    u = NaN;
end

end