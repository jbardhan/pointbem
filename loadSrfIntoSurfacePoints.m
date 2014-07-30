function surfdata = loadSrfIntoSurfacePoints(srf)

[meshBase,rootDir] = readsrf(srf);
meshData = readmesh(meshBase,1);
[panelCentroids,panelNormals,panelAreas] = genmeshcolloc(meshData);

vertexWeights = zeros(size(meshData.vert,1),1);
for i=1:size(meshData.face,1)
  vertexWeights(meshData.face(i,1:3)) = vertexWeights(meshData.face(i,1:3)) + panelAreas(i)/3.0;
end

surfdata = struct('points',meshData.vert,'weights',vertexWeights, ...
		  'anglecoords', 0,'normals',meshData.normals);

% to plot: 
% v = meshData.vert;
% n = meshData.normals;
% quiver3(v(:,1),v(:,2),v(:,3),n(:,1),n(:,2),n(:,3));
