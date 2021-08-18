function [optimXArray, optimCostArray, USample, optimTime] = optimizeBpf_fmincon( targetBpf, targetCov, initX, evalLimitFunVal_singleFun, nLimitFun, nV, vStdArray, costFun, constFun, xLowerBound, xUpperBound, solutionStepRatio, activeAndFailRatio, solutionChangeRatioTol, optimOption, isViolationPenalty )

import function.*

nSample = ceil( (1-targetBpf)/targetBpf/targetCov^2 );
nSampleFail = ceil( nSample * targetBpf );
nSampleActive = ceil( nSampleFail * activeAndFailRatio );
USample = normrnd( 0, 1, [nSample, nV]);

optimUpperBound = [xUpperBound inf( 1,nLimitFun ) inf(1,nSampleActive*nLimitFun)]; 
optimLowerBound = [xLowerBound -inf( 1,nLimitFun ) zeros(1,nSampleActive*nLimitFun)]; 

nX = length( xLowerBound );


if nargin < 16
    isViolationPenalty = false; % default
end

if ~isViolationPenalty
    nDecisionVariable = nX + nLimitFun + nLimitFun * nSampleActive;
    Aineq = zeros( nLimitFun, nDecisionVariable);
    for iLimitFunIndex = 1:nLimitFun
        iAineq = zeros( 1, nDecisionVariable );
        iAineq( nX + iLimitFunIndex ) = 1;
        iAineq( nX + nLimitFun + nSampleActive*(iLimitFunIndex-1)+(1:nSampleActive) ) = 1/targetBpf/nSample;

        Aineq( iLimitFunIndex, : ) = iAineq;
    end
    bineq = zeros( nLimitFun, 1 );
else 
    Aineq = []; bineq = [];
    costFun = @(x) costFun( x, nSample );
end


optimXArray = zeros( nX, 0 ); optimCostArray = zeros( 1, 0 );
nInteration = 0;
iSolutionChangeRatio = 100;
iX = initX;
tic
while iSolutionChangeRatio > solutionChangeRatioTol
    nInteration = nInteration + 1;
    iXOld = iX;
    iV = repmat( vStdArray(:).', nSample, 1 ) .* USample + repmat( iXOld(:).', nSample, 1 );
    
    iLimitFunValSortActive = zeros( nSampleActive, nLimitFun );
    iLimitFunValSortActiveIndex = zeros( nSampleActive, nLimitFun );
    iInitialSolution = zeros( nX + nLimitFun + nLimitFun*nSampleActive, 1 ); iInitialSolution(1:nX) = iXOld;
    for ijLimitFunIndex = 1:nLimitFun
        ijLimitFunVal = evalLimitFunVal_singleFun( iV, false, ijLimitFunIndex );
        [ijLimitFunValSort, ijLimitFunValSortInd] = sort( ijLimitFunVal, 'descend' );

        ijLimitFunValSortActive = ijLimitFunValSort( 1:nSampleActive );
        iLimitFunValSortActive( :, ijLimitFunIndex ) = ijLimitFunValSortActive;
        
        ijLimitFunValSortActiveIndex = ijLimitFunValSortInd( 1:nSampleActive );
        iLimitFunValSortActiveIndex( :, ijLimitFunIndex ) = ijLimitFunValSortActiveIndex;

        ijZ0andZInit = evalZValues( ijLimitFunValSortActive, nSample, targetBpf );
        iInitialSolution( nX + ijLimitFunIndex ) = ijZ0andZInit(1);
        iInitialSolution( nX + nLimitFun + nSampleActive * ( ijLimitFunIndex-1 ) + (1:nSampleActive ) ) = ijZ0andZInit(2:end);
    end
    

    
    [iXZ,iFval] = fmincon( @(x) costFun( x ), iInitialSolution, Aineq, bineq, [],[],optimLowerBound,optimUpperBound,@(x) constFun( x, USample, iLimitFunValSortActiveIndex ), optimOption);

    iX = iXZ( 1:nX );
    iX = iX(:)';
    iSolutionStepRatio = solutionStepRatio^( nInteration-1 );
    iX = iXOld + iSolutionStepRatio * (iX - iXOld);
    iSolutionChangeRatio = norm( iX - iXOld ) / norm( iX );
    
    optimXArray = [optimXArray iX(:)];
    optimCostArray = [optimCostArray iFval];
end
optimTime = toc;