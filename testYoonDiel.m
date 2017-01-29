global kappa;

radius   = 2.5;
origin   = [0 0 0];
epsIn = 1;
epsOut = 80;
conv_factor = 332.112;  % so our units (energy) come out in kilocalories per mole
pqr = struct('xyz',origin,'q',1,'R',1);
dG_exact = 0.5 * conv_factor * (1/epsOut - 1/epsIn) * 1/radius;

density = 2; % vertices per Angstrom squared
numpoints(i) = floor(density * pi * radius^2);
sph = makeSphereSurface(origin, radius, numpoints(i));
yoon = makeBemYoonDielMatrices(sph, pqr, epsIn, epsOut);
dG_solv = 0.5 * conv_factor * yoon.C * (yoon.A\yoon.B);
