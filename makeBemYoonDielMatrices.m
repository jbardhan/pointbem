function bem = makeBemYoonDielMatrices(dielSrf, pqr, epsIn, epsOut)

npDiel = length(dielSrf.weights);
nCharges = length(pqr.q);

dielChargeOp = makeSurfaceToChargeOperators(dielSrf, pqr);

Idiel = eye(npDiel);

dielDielOp = makeSurfaceToSurfaceLaplaceOperators(dielSrf);

			     
A = [(Idiel/2 + dielDielOp.K) -dielDielOp.V ;
     (Idiel/2 - dielDielOp.K) ((epsIn/epsOut)*dielDielOp.V)];

B = [dielChargeOp.phiCoul/epsIn;
     zeros(npDiel,nCharges);];

C = 4*pi*[-dielChargeOp.dlpToCharges dielChargeOp.slpToCharges];

bem = struct('B', B, 'A', A, 'C', C,...
	     'dielChargeOp',dielChargeOp,...
	     'dielDielOp',dielDielOp);