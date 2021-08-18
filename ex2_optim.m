%{
11 May 21; 12 May 21
Ji-Eun Byun

Ex2 - optimization
%}

clear; close all;
import ex2.*
import function.*
rng(2)

%% Problem
nX = 4;
xLowerBound = [3.175 15 200 3.175]; xUpperBound = [10 254 220 10];

nLimitFun = 5;
nV = 4;
vStdArray = [0.1693 0.1693 0.0107 0.0107];

evalLimitFunVal = @(vRowArray) limitFunEx2( vRowArray );
evalLimitFunVal_singleFun = @(vRowArray, gradientTrueOrFalse, limitFunIndex) limitFunEx2_singleFun( vRowArray, gradientTrueOrFalse, limitFunIndex );

%% Optimization setting
solutionStepRatio = 0.9;
activeAndFailRatio = 1.2;
solutionChangeRatioTol = 1e-2;
optimOption = optimoptions('fmincon', 'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, 'MaxFunctionEvaluations', 1e5, 'MaxIterations', 1e4 );

%% Problem-specific functions
funEvalCost = @(x) evalCostEx2( x );
constFun = @( x, uRowArray, activeSampleIndexArray ) constraintAndGradient( x, uRowArray, activeSampleIndexArray, evalLimitFunVal_singleFun, nX, vStdArray );

%% Reference solution
nManySample = 1e6;
manyUSample = randn( nManySample,nV );

optimRoundName = 'Reference';   
referenceX = [5.72 2e2 2.11e2 6.25].';
[referencePf, referenceBpf] = summarizeOptimResult( referenceX, manyUSample, [], funEvalCost, evalLimitFunVal, vStdArray, [], [], optimRoundName );
referenceCost = funEvalCost( referenceX );

targetBpf = max( referenceBpf );


%% Find feasible solution by Lagrangian relaxation
targetCov0 = 0.2;
initX = 0.5 * ( xLowerBound + xUpperBound );
penalty = 10;
costFun_Lagrangian = @(x, nSample ) costAndGradientEx2( x, nLimitFun, targetBpf, penalty, nSample );

[optimXArray0, optimCostArray0, USample0, optimTime0] = optimizeBpf_fmincon( targetBpf, targetCov0, initX, evalLimitFunVal_singleFun, nLimitFun, ...
                                                                             nV, vStdArray, costFun_Lagrangian, constFun, xLowerBound, xUpperBound, ...
                                                                             solutionStepRatio, activeAndFailRatio, solutionChangeRatioTol, optimOption, true );
                                                                                                                                                                                                                           
%% Optimization with low c.o.v.
targetCov1 = 0.2;
initX = optimXArray0(:,end).';
costFun = @(x) costAndGradientEx2( x );

[optimXArray1, optimCostArray1, USample1, optimTime1] = optimizeBpf_fmincon( targetBpf, targetCov1, initX, evalLimitFunVal_singleFun, nLimitFun, ...
                                                                             nV, vStdArray, costFun, constFun, xLowerBound, xUpperBound, ...
                                                                             solutionStepRatio, activeAndFailRatio, solutionChangeRatioTol, optimOption );
                                                                                                                    
%% Optim with high c.o.v.
% targetCov2 = 0.1;
targetCov2 = 0.05;
initX = optimXArray1(:,end).';

[optimXArray2, optimCostArray2, USample2, optimTime2] = optimizeBpf_fmincon( targetBpf, targetCov2, initX, evalLimitFunVal_singleFun, nLimitFun, ...
                                                                             nV, vStdArray, costFun, constFun, xLowerBound, xUpperBound, ...
                                                                             solutionStepRatio, activeAndFailRatio, solutionChangeRatioTol, optimOption );
    

%% Result
evalLimitFunVal = @(vRowArray) limitFunEx2( vRowArray );
optimRoundName = 'Penalty function';                                                                                                                   
[optimPf0, optimBpf0, limitFunValArray0, z0_0, zArray0, z0Index0] = summarizeOptimResult( optimXArray0, manyUSample, optimTime0, @(x) ex2.evalCostEx2( x ), @(vRowArray) ex2.limitFunEx2( vRowArray ), vStdArray, max(referencePf), targetBpf, optimRoundName ); 

optimRoundName = 'Low c.o.v.';                                                                                                                   
[optimPf1, optimBpf1, limitFunValArray1, z0_1, zArray1, z0Index1] = summarizeOptimResult( optimXArray1, manyUSample, optimTime1, @(x) ex2.evalCostEx2( x ), @(vRowArray) ex2.limitFunEx2( vRowArray ), vStdArray, max(referencePf), targetBpf, optimRoundName );

optimRoundName = 'High c.o.v.';                                                                                                                      
[optimPf2, optimBpf2, limitFunValArray2, z0_2, zArray2, z0Index2, vArray2] = summarizeOptimResult( optimXArray2, manyUSample, optimTime2,  @(x) ex2.evalCostEx2( x ), @(vRowArray) ex2.limitFunEx2( vRowArray ), vStdArray, max(referencePf), targetBpf, optimRoundName );


save ex2_optim