function bem = makeBemSternMatrices(dielSrf, sternSrf, pqr, epsIn, ...
				    epsOut, kappaArg)
global kappa;
kappa = kappaArg;

npDiel = length(dielSrf.weights)
npStern = length(sternSrf.weights)
nCharges = length(pqr.q);

dielChargeOp = makeSurfaceToChargeOperators(dielSrf, pqr);

Idiel = eye(npDiel);
Istern = eye(npStern);

% src = diel, target = diel 
dielDielOp = makeSurfaceToSurfaceLaplaceOperators(dielSrf);

% src = diel, target = stern
dielSternOp = makeSurfaceToSurfaceLaplaceOperators(dielSrf, ...
						  sternSrf);

% src = stern, target = diel
sternDielOp = makeSurfaceToSurfaceLaplaceOperators(sternSrf, ...
						  dielSrf);

% LAPLACE kernel for src = stern, target = stern
sternSternLaplace = makeSurfaceToSurfaceLaplaceOperators(sternSrf);

% YUKAWA kernel for src = stern, target = stern
sternSternYukawa = makeSurfaceToSurfaceYukawaOperators(sternSrf);

Z1 = zeros(npDiel, npStern);
Z2 = zeros(npStern, npDiel); % could just use Z1' but this is
                             % clearer
			     
A = [(Idiel/2 + dielDielOp.K) -dielDielOp.V Z1 Z1;
     (Idiel/2 - dielDielOp.K) ((epsIn/epsOut)*dielDielOp.V) ...
     sternDielOp.K -sternDielOp.V;
     -dielSternOp.K ((epsIn/epsOut)*dielSternOp.V)...
     (Istern/2+sternSternLaplace.K) -sternSternLaplace.V;
     Z2 Z2 (Istern/2-sternSternYukawa.K) ((epsOut/epsOut)* ...
					  sternSternYukawa.V)];


B = [dielChargeOp.phiCoul/epsIn;
     zeros(npDiel,nCharges);
     zeros(npStern,nCharges);
     zeros(npStern,nCharges);];

C = [-dielChargeOp.dlpToCharges dielChargeOp.slpToCharges ...
     zeros(nCharges,npStern) zeros(nCharges, npStern)];

bem = struct('B', B, 'A', A, 'C', C,...
	     'dielChargeOp',dielChargeOp,...
	     'dielDielOp',dielDielOp,...
	     'dielSternOp',dielSternOp,...
	     'sternDielOp',sternDielOp,...
	     'sternSternLaplace',sternSternLaplace,...
	     'sternSternYukawa',sternSternYukawa);
