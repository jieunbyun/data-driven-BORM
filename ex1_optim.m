%{
21 Apr 21; 12 May 21
Ji-Eun Byun

Ex1 - optimization
%}

clear; close all;
import ex1.*
import function.*
rng(100)

%% Problem
nX = 2;
xLowerBound = [0 0]; xUpperBound = [3.7 4.0];

nLimitFun = 2;
nV = 2;
vStdArray = 0.1 * ones( 1, nV );

evalLimitFunVal = @(vRowArray) limitFunEx1( vRowArray );
evalLimitFunVal_singleFun = @(vRowArray, gradientTrueOrFalse, limitFunIndex) limitFunEx1_singleFun( vRowArray, gradientTrueOrFalse, limitFunIndex );

%% Optimization setting
solutionStepRatio = 0.9;
activeAndFailRatio = 1.2;
solutionChangeRatioTol = 1e-2;
optimOption = optimoptions('fmincon', 'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);

%% Problem-specific functions
funEvalCost = @(x) evalCostEx1( x );
constFun = @(x, uRowArray, activeSampleIndexArray) constraintAndGradient( x, uRowArray, activeSampleIndexArray, evalLimitFunVal_singleFun, nX, vStdArray );


%% Reference solution
nManySample = 1e6;
manyUSample = randn( nManySample,nV );

optimRoundName = 'Reference';   
referenceX = [2.81 3.28].';
[referencePf, referenceBpf] = summarizeOptimResult( referenceX, manyUSample, [], funEvalCost, evalLimitFunVal, vStdArray, [], [], optimRoundName );

targetBpf = max(referenceBpf);

%% Find feasible solution by Lagrangian relaxation
targetCov0 = 0.2;
initX = [2 2];
penalty = 10;
costFun_Lagrangian = @(x,nSample) costAndGradientEx1( x, nLimitFun, targetBpf, penalty, nSample );

[optimXArray0, optimCostArray0, USample0, optimTime0] = optimizeBpf_fmincon( targetBpf, targetCov0, initX, evalLimitFunVal_singleFun, nLimitFun, ...
                                                                             nV, vStdArray, costFun_Lagrangian, constFun, xLowerBound, xUpperBound, ...
                                                                             solutionStepRatio, activeAndFailRatio, solutionChangeRatioTol, optimOption, true );
                                                                                                                                                                                                                           
%% Optimization with low c.o.v.
targetCov1 = 0.2;
initX = optimXArray0(:,end).';
costFun = @(x) costAndGradientEx1( x );

[optimXArray1, optimCostArray1, USample1, optimTime1] = optimizeBpf_fmincon( targetBpf, targetCov1, initX, evalLimitFunVal_singleFun, nLimitFun, ...
                                                                             nV, vStdArray, costFun, constFun, xLowerBound, xUpperBound, ...
                                                                             solutionStepRatio, activeAndFailRatio, solutionChangeRatioTol, optimOption );
                                                                                                                    
%% Optim with high c.o.v.
targetCov2 = 0.05;
initX = optimXArray1(:,end).';

[optimXArray2, optimCostArray2, USample2, optimTime2] = optimizeBpf_fmincon( targetBpf, targetCov2, initX, evalLimitFunVal_singleFun, nLimitFun, ...
                                                                             nV, vStdArray, costFun, constFun, xLowerBound, xUpperBound, ...
                                                                             solutionStepRatio, activeAndFailRatio, solutionChangeRatioTol, optimOption );
                                                                                                                    
%% Result
optimRoundName = 'Penalty function';                                                                                                                   
[optimPf0, optimBpf0, limitFunValArray0, z0_0, zArray0, z0Index0] = summarizeOptimResult( optimXArray0, manyUSample, optimTime0, funEvalCost, evalLimitFunVal, vStdArray, referencePf, targetBpf, optimRoundName ); 

optimRoundName = 'Low c.o.v.';                                                                                                                   
[optimPf1, optimBpf1, limitFunValArray1, z0_1, zArray1, z0Index1] = summarizeOptimResult( optimXArray1, manyUSample, optimTime1, funEvalCost, evalLimitFunVal, vStdArray, referencePf, targetBpf, optimRoundName );

optimRoundName = 'High c.o.v.';                                                                                                             
[optimPf2, optimBpf2, limitFunValArray2, z0_2, zArray2, z0Index2, vArray2] = summarizeOptimResult( optimXArray2, manyUSample, optimTime2, funEvalCost, evalLimitFunVal, vStdArray, referencePf, targetBpf, optimRoundName );


% Plot failure samples
[~, ~, ~, ~, zArray2_failPlot, ~, vArray2_failPlot] = summarizeOptimResult( optimXArray2, USample2, optimTime2, funEvalCost, evalLimitFunVal, vStdArray, referencePf, targetBpf, optimRoundName );

markerSize = 7.5; fontsize_tick = 12; fontName = 'times new roman'; fontsize_label = 16;
figure;
plot( vArray2_failPlot(zArray2_failPlot>0,1), vArray2_failPlot(zArray2_failPlot>0,2), '.', 'markersize', markerSize )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
xlabel( '{\it V}_1', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it V}_2', 'fontname', fontName, 'fontsize', fontsize_label )

saveas( gcf, 'figure/ex1_failedSamples.emf' )
saveas( gcf, 'figure/ex1_failedSamples.pdf' )


save ex1_optim