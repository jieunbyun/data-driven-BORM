function [optimXArray, optimCostArray, USample, optimTime] = optimizeBpf_fmincon_activeLimitStateFun( targetBpf, targetCov, initX, evalLimitFunVal_singleFun, nLimitFun, nV, vStdArray, costFun, constFun, xLowerBound, xUpperBound, solutionStepRatio, activeAndFailRatio, solutionChangeRatioTol, limitStateActiveRatioWithBpf, optimOption, isViolationPenalty )

import function.*

nSample = ceil( (1-targetBpf)/targetBpf/targetCov^2 );
nSampleFail = ceil( nSample * targetBpf );
nSampleActive = ceil( nSampleFail * activeAndFailRatio );
USample = normrnd( 0, 1, [nSample, nV]);

nX = length( xLowerBound );


if nargin < 17
    isViolationPenalty = false; % default
end

if isViolationPenalty
    iAineq = []; ibineq = [];
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
    iZ0andZInit = zeros( 1+nSampleActive, nLimitFun );
    iLimitFunBpf = zeros( 1, nLimitFun );
    for ijLimitFunIndex = 1:nLimitFun
        
        ijLimitFunVal = evalLimitFunVal_singleFun( iV, false, ijLimitFunIndex );
        iLimitFunBpf( ijLimitFunIndex ) = evalBpf( ijLimitFunVal );
        
        [ijLimitFunValSort, ijLimitFunValSortInd] = sort( ijLimitFunVal, 'descend' );

        ijLimitFunValSortActive = ijLimitFunValSort( 1:nSampleActive );
        iLimitFunValSortActive(:, ijLimitFunIndex ) = ijLimitFunValSortActive;

        ijLimitFunValSortActiveIndex = ijLimitFunValSortInd( 1:nSampleActive );
        iLimitFunValSortActiveIndex(:, ijLimitFunIndex ) = ijLimitFunValSortActiveIndex;

        ijZ0andZInit = evalZValues( ijLimitFunValSortActive, nSample, targetBpf );
        iZ0andZInit(:, ijLimitFunIndex) = ijZ0andZInit;
        
    end
    
    iActiveLimitFunIndex = find( iLimitFunBpf/targetBpf > limitStateActiveRatioWithBpf );
    if isempty( iActiveLimitFunIndex )
        [~,iActiveLimitFunIndex] = max( iLimitFunBpf );
    end
    iNActiveLimitFun = length( iActiveLimitFunIndex );
    iLimitFunValSortActiveIndex = iLimitFunValSortActiveIndex( :,iActiveLimitFunIndex );
    
    iInitialSolution = zeros( nX + iNActiveLimitFun + iNActiveLimitFun*nSampleActive, 1 ); iInitialSolution(1:nX) = iXOld;
    iZ0andZInit = iZ0andZInit(:,iActiveLimitFunIndex);
    for ijActiveLimitFunIndex = 1:iNActiveLimitFun
        ijZ0andZInit = iZ0andZInit( :, ijActiveLimitFunIndex );
        iInitialSolution( nX + ijActiveLimitFunIndex ) = ijZ0andZInit(1);
        iInitialSolution( nX + iNActiveLimitFun + nSampleActive * ( ijActiveLimitFunIndex-1 ) + (1:nSampleActive ) ) = ijZ0andZInit(2:end);
    end           
    
    
    if ~isViolationPenalty
        iNDecisionVariable = nX + iNActiveLimitFun + iNActiveLimitFun * nSampleActive;
        iAineq = zeros( iNActiveLimitFun, iNDecisionVariable);
        for ijActiveLimitFunIndex = 1:iNActiveLimitFun
            ijAineq = zeros( 1, iNDecisionVariable );
            ijAineq( nX + ijActiveLimitFunIndex ) = 1;
            ijAineq( nX + iNActiveLimitFun + nSampleActive*(ijActiveLimitFunIndex-1)+(1:nSampleActive) ) = 1/targetBpf/nSample;

            iAineq( ijActiveLimitFunIndex, : ) = ijAineq;
        end
        ibineq = zeros( iNActiveLimitFun, 1 );
        
        iCostFun = costFun;
    else
        iCostFun = @(x) costFun( x,iNActiveLimitFun, nSample );
    end
    
    iOptimUpperBound = [xUpperBound inf( 1,iNActiveLimitFun ) inf(1,nSampleActive*iNActiveLimitFun)]; 
    iOptimLowerBound = [xLowerBound -inf( 1,iNActiveLimitFun ) zeros(1,nSampleActive*iNActiveLimitFun)]; 
    
    
    [iXZ,iFval] = fmincon( @(x) iCostFun( x ), iInitialSolution, iAineq, ibineq, [],[],iOptimLowerBound,iOptimUpperBound,@(x) constFun( x, USample, iLimitFunValSortActiveIndex, iActiveLimitFunIndex ), optimOption);

    iX = iXZ( 1:nX );
    iX = iX(:)';
    iSolutionStepRatio = solutionStepRatio^( nInteration-1 );
    iX = iXOld + iSolutionStepRatio * (iX - iXOld);
    iSolutionChangeRatio = norm( iX - iXOld ) / norm( iX );
    
    optimXArray = [optimXArray iX(:)];
    optimCostArray = [optimCostArray iFval];
end
optimTime = toc;