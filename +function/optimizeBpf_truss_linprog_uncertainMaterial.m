function [optimXArray, optimCostArray, loadSampleArraySelected, optimTime, yieldStrengthSampleArray] = optimizeBpf_truss_linprog_uncertainMaterial( targetBpf, targetCov, activeAndFailRatio, solutionStepRatio, solutionChangeRatioTol, loadSampleArray, xLowerBound, xUpperBound, xId2MemberId, memberForceByUnitWeight, memberLengthArray, steelUnitWeight_knM3, steelYieldStrength_knM2, steelYieldStrengthCov )

import ex3.*

nSample = ceil( (1-targetBpf)/targetBpf/targetCov^2 );
loadSampleArraySelected = loadSampleArray( :, 1:nSample );
nSampleFail = ceil( nSample * targetBpf );
nSampleActive = ceil( nSampleFail * activeAndFailRatio );

nX = length( xId2MemberId );
nMember = size( loadSampleArray, 1 );

yieldStrengthSampleArray = steelYieldStrength_knM2 + steelYieldStrength_knM2 * steelYieldStrengthCov *normrnd( 0, 1, [nMember, nSample] );

nZ0 = nMember; nZ = nMember*nSampleActive; nS = nMember*nSampleActive;
lowerBound = [xLowerBound -inf( 1,nZ0 ) zeros( 1,nZ ) -inf( 1,nS )]; upperBound = [xUpperBound inf( 1,nZ0 ) inf( 1,nZ ) inf( 1,nS )]; 

costFun = getCostFunction( nSampleActive, xId2MemberId, memberLengthArray, steelUnitWeight_knM3 );
[Aineq_z0, bineq_z0] = getInequalityConstraint_z0( nX, nMember, nSample, nSampleActive, targetBpf );

optimXArray = zeros( nX, 0 ); optimCostArray = zeros( 1, 0 );
nInteration = 0;
iSolutionChangeRatio = 100;
iX = xLowerBound;
tic
while iSolutionChangeRatio > solutionChangeRatioTol
    nInteration = nInteration + 1;
    iXOld = iX;
    
    iLimitStateFunValArray = evalLimitFun_uncertainMaterial( iXOld, xId2MemberId, loadSampleArraySelected, memberForceByUnitWeight, memberLengthArray, steelUnitWeight_knM3, yieldStrengthSampleArray );
    [~, iActiveSampleArrayIndex] = findActiveSample( iLimitStateFunValArray, nSampleActive );
    [iAineq_z, iBineq_z] = getInequalityConstraint_z_uncertainMaterial( xId2MemberId, yieldStrengthSampleArray, iActiveSampleArrayIndex );
    [iAineq_s, iBineq_s] = getInequalityConstraint_s( xId2MemberId, loadSampleArraySelected, iActiveSampleArrayIndex, memberForceByUnitWeight, memberLengthArray, steelUnitWeight_knM3 );
    
    [iXZ, iFval] = linprog( costFun, [Aineq_z0; iAineq_z; iAineq_s], [bineq_z0; iBineq_z; iBineq_s], [], [], lowerBound, upperBound );
    
    iX = iXZ( 1:nX );
    iX = iX(:)';
    iSolutionStepRatio = solutionStepRatio^( nInteration-1 );
    iX = iXOld + iSolutionStepRatio * (iX - iXOld);
    iSolutionChangeRatio = norm( iX - iXOld ) / norm( iX );
    
    optimXArray = [optimXArray iX(:)];
    optimCostArray = [optimCostArray iFval];
end
optimTime = toc;






