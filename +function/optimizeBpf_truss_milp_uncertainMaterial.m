function [optimXArray, optimCostArray, loadSampleArraySelected, optimTime] = optimizeBpf_truss_milp_uncertainMaterial( targetBpf, targetCov, activeAndFailRatio, solutionChangeRatioTol, loadSampleArray, xSolutionSet, xId2MemberId, memberForceByUnitWeight, memberLengthArray, steelUnitWeight_knM3, steelYieldStrength_knM2, steelYieldStrengthCov )

import ex3.*

nSample = ceil( (1-targetBpf)/targetBpf/targetCov^2 );
loadSampleArraySelected = loadSampleArray( :, 1:nSample );
nSampleFail = ceil( nSample * targetBpf );
nSampleActive = ceil( nSampleFail * activeAndFailRatio );

nX = length( xId2MemberId );
nMember = size( loadSampleArray, 1 );
nSolution = length( xSolutionSet );

yieldStrengthSampleArray = steelYieldStrength_knM2 + steelYieldStrength_knM2 * steelYieldStrengthCov *normrnd( 0, 1, [nMember, nSample] );

nZ0 = nMember; nZ = nMember*nSampleActive; nS = nMember*nSampleActive;
lowerBound = [zeros( 1,nX*nSolution ) -inf( 1,nZ0 ) zeros( 1,nZ ) -inf( 1,nS )]; upperBound = [ones( 1,nX*nSolution ) inf( 1,nZ0 ) inf( 1,nZ ) inf( 1,nS )]; 

costFun = getCostFunction_milp( nSampleActive, xId2MemberId, memberLengthArray, steelUnitWeight_knM3, xSolutionSet );
[Aineq_z0, bineq_z0] = getInequalityConstraint_z0( nX*nSolution, nMember, nSample, nSampleActive, targetBpf );
[Aeq_x, beq_x] = getEqualityConstraint_x_milp( nX, nSolution, nZ0, nZ, nS );
integerDecisionVariableIndex = 1:(nX*nSolution);

optimXArray = zeros( nX, 0 ); optimCostArray = zeros( 1, 0 );
nInteration = 0;
iSolutionChangeRatio = 100;
iX = min( xSolutionSet ) * ones( 1, nX );
tic
while iSolutionChangeRatio > solutionChangeRatioTol
    nInteration = nInteration + 1;
    iXOld = iX;
    
    iLimitStateFunValArray = evalLimitFun_uncertainMaterial( iXOld, xId2MemberId, loadSampleArraySelected, memberForceByUnitWeight, memberLengthArray, steelUnitWeight_knM3, yieldStrengthSampleArray );
    [~, iActiveSampleArrayIndex] = findActiveSample( iLimitStateFunValArray, nSampleActive );
    [iAineq_z, iBineq_z] = getInequalityConstraint_z_uncertainMaterial_milp( xId2MemberId, yieldStrengthSampleArray, iActiveSampleArrayIndex, xSolutionSet );
    [iAineq_s, iBineq_s] = getInequalityConstraint_s_milp( xId2MemberId, loadSampleArraySelected, iActiveSampleArrayIndex, memberForceByUnitWeight, memberLengthArray, steelUnitWeight_knM3, xSolutionSet );
    
    [iXZ, iFval] = intlinprog( costFun, integerDecisionVariableIndex, [Aineq_z0; iAineq_z; iAineq_s], [bineq_z0; iBineq_z; iBineq_s], Aeq_x, beq_x, lowerBound, upperBound );
    
    iX = zeros( 1, nX );
    for jXIndex = 1:nX
        jXLocationArray = nSolution * (jXIndex-1) + (1:nSolution);
        jX = iXZ( jXLocationArray );
        jSolution = xSolutionSet( jX > 0+1e-3 );
        iX( jXIndex ) = jSolution;
    end
    
    iSolutionChangeRatio = norm( iX - iXOld ) / norm( iX );
    
    optimXArray = [optimXArray iX(:)];
    optimCostArray = [optimCostArray iFval];
end
optimTime = toc;






