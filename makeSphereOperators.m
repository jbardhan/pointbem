function operstruct = makeSphereOperators(origin, radius, numpoints)

[pnts,wts,angs,normals] = getSphPoints(origin, radius, numpoints);
[V,K] = genPointMatrices(pnts,normals,wts);

if 0
  [V_eigenvecs,V_eigenvals]= eig(V);
  lambda=diag(V_eigenvals);
  [j1,I] = sort(real(lambda),'descend');
  V_eigenvals = lambda(I);
  V_eigenvecs = V_eigenvecs(:,I);
  
  [K_eigenvecs,K_eigenvals]= eig(K);
  lambda=diag(K_eigenvals);
  [j1,I] = sort(real(lambda),'ascend');
  K_eigenvals = lambda(I);
  K_eigenvecs = K_eigenvecs(:,I);
else
  V_eigenvecs = 0; V_eigenvals = 0;
  K_eigenvecs = 0; K_eigenvals = 0;
end

operstruct = struct('points',pnts,'weights',wts, 'anglecoords', ...
		    angs,'normals',normals, 'V', V, 'K', K, ...
		    'V_eigenvals',V_eigenvals, 'V_eigenvecs',V_eigenvecs, ...
		    'K_eigenvals',K_eigenvals, 'K_eigenvecs',K_eigenvecs);