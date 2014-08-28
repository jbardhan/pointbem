function [V, K, W] = genPointYukawaMatrices(points, normals, weights,...
					 destpoints, destnormals, destweights)
global kappa;
printInfo = 1;

if nargin < 6
  destpoints = points;
  destnormals = normals;
  destweights = weights;
end

Nsrc = length(weights);
Ndest = length(destweights);
V = zeros(Ndest, Nsrc);
K = V;
W = V;
for i=1:Ndest
  for j=1:Nsrc
	 vecr = destpoints(i,:)-points(j,:);
	 r = norm(vecr);
	 if r > 1e-6
		V(i,j) = weights(j)* exp(-kappa*r)/4/pi/r;
		K(i,j) = weights(j)* vecr*normals(j,:)'*(1+kappa*r)*exp(-kappa*r)/4/pi/r/r/r;

%		n0 = destnormals(i,:);
%		n  = normals(j,:);
%		cosTheta = vecr*n'/r;
%		cosTheta0 = vecr*n0'/r;
%		d2Fdn0dn = weights(j)*((n0*n')-3*cosTheta*cosTheta0)/4/pi/r/r/r;
%		W(i,j) = weights(j)*((1+kappa*r)*exp(-kappa*r)*d2Fdn0dn ...
%			 - kappa^2 * exp(-kappa*r)*cosTheta*cosTheta0/4/pi/r);
	 end
  end
end

if printInfo 
  fprintf('genPointYukawaMatrices: calculating hypersingular');
  fprintf(' operator using Calderon formula \n\tW = V^{-1} (1/4 * (I - K^2))\n');
end

[V_eigenvecs,V_eigenvals]= eig(V);
lambda=diag(V_eigenvals);
[j1,I] = sort(real(lambda),'descend');
V_eigenvals = lambda(I);
V_eigenvecs = V_eigenvecs(:,I);
W = V_eigenvecs * (diag(1./V_eigenvals) * (V_eigenvecs'*(.25*(eye(size(K,1))-K^2))));
