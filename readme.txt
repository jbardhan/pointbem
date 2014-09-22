test scripts
* testConvergenceYoonDiel.m
* testHypersingular.m
* testSrfOnSurfacePoints.m
* testStern.m
* testpoint.m

* makeComparison.m

* geometry/ directory
* loadConstants.m

------------------------------------
------------------------------------
fileIO

* readmesh.m
* readpqr.m
* readsrf.m
* parsePQRentry.m

------------------------------------
------------------------------------
constructing test problems

* addPqrGridSpherePlane.m

------------------------------------
------------------------------------
simple primitives and exact results

* besselBackwardRecurrence.m
* besselBackwardRecurrenceDeriv.m
* besselCai.m
* besselDerivative.m
* besselForwardRecurrence.m
* besselForwardRecurrenceDeriv.m

* computeLaplaceEigenvalues.m
* computeLaplaceEigenvaluesConc.m
* computeYukawaEigenvalues.m

* convertToSph.m
* getSphPoints.m

* computeBnm_exact.m
* computeEnm.m
* computePot.m
* doAnalytical.m

------------------------------------
------------------------------------
pointBEM

* loadSrfIntoSurfacePoints.m
* getPointCoulomb.m
* genPointLaplaceMatrices.m
* genPointYukawaMatrices.m

* makeBemEcfQualMatrices.m
* makeBemSternMatrices.m
* makeBemYoonDielMatrices.m
* makeBemYoonLPBMatrices.m

* makeSurfaceToChargeOperators.m
* makeSurfaceToSurfaceLaplaceOperators.m
* makeSurfaceToSurfaceOperators.m
* makeSurfaceToSurfaceYukawaOperators.m

------------------------------------
------------------------------------
general (panel) BEM

* calcp.m
* genmeshcolloc.m

------------------------------------
------------------------------------
Special case sphere problem:

* makeSphereChargeDistribution.m
* makeSphereSurface.m
* makeSphereLaplaceOperators.m
* makeSphereOperators.m
* makeSphereYukawaOperators.m

