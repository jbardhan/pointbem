function bem = makeBemYoonLPBMatrices(dielSrf, pqr, epsIn, epsOut, kappaArg)
global kappa;
kappa = kappaArg;

npDiel = length(dielSrf.weights);
nCharges = length(pqr.q);

dielChargeOp = makeSurfaceToChargeOperators(dielSrf, pqr);

Idiel = eye(npDiel);

dielDielLaplace = makeSurfaceToSurfaceLaplaceOperators(dielSrf);
dielDielYukawa =  makeSurfaceToSurfaceYukawaOperators(dielSrf);
			     
A = [(Idiel/2 + dielDielLaplace.K) -dielDielLaplace.V ;
     (Idiel/2 - dielDielYukawa.K) ((epsIn/epsOut)*dielDielYukawa.V)];

B = [dielChargeOp.phiCoul/epsIn;
     zeros(npDiel,nCharges);];

C = 4*pi*[-dielChargeOp.dlpToCharges dielChargeOp.slpToCharges];

bem = struct('B', B, 'A', A, 'C', C,...
	     'dielChargeOp',dielChargeOp,...
	     'dielDielOp',dielDielOp);