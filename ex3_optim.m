%{
26 Apr 21
Ji-Eun Byun

Ex3 - optimization
%}

clear;
import ex3.*
import function.*
rng(1)

%% Truss structure
load data/truss50m_member10m.mat
nX = 6;
xId2MemberId = cell( nX, 1 );
xId2MemberId{1} = [1 nMember-1];
xId2MemberId{2} = [2 nMember];
xId2MemberId{3} = 3:5:(nMember-1);
xId2MemberId{4} = 4:5:(nMember-2);
xId2MemberId{5} = sort( [5:5:(nMember-1) 6:5:(nMember-1)] );
xId2MemberId{6} = 7:5:(nMember-1);

memberForceByUnitWeight = evalMemberForceByUnitWeight( Connectivity, memberForceArray );
memberLength_m = evalMemberLength( Connectivity, Coord );

xLowerBound_m2 = 1.0e3 * 1e-6 * ones(1,nX);
xUpperBound_m2 = 4.0e3 * 1e-6 * ones(1,nX);
steelUnitWeight_kgM3 = 7950;
g_mS2 = 9.8; 
steelUnitWeight_knM3 = steelUnitWeight_kgM3 * g_mS2 * 1e-3;
steelYieldStrength_knM2 = 250*1e3;

%% Problem
targetPf = 1e-2;
tau = getBufferedTailIndexFromPf( targetPf );
% targetBpf = targetPf * 2.65;
targetBpf = targetPf * tau;
nLimitFun = nMember;
nV = nMember;

load data/loadData.mat
solutionStepRatio = 0.9;
activeAndFailRatio = 1.2;
solutionChangeRatioTol = 1e-2;

targetCov = 0.05;
%% Optimization
% LP
[optimXArray_m2, optimCostArray, loadSampleArraySelected, optimTime] = optimizeBpf_truss_linprog( targetBpf, targetCov, activeAndFailRatio, solutionStepRatio, solutionChangeRatioTol, maxMemberForceInKNArray, xLowerBound_m2, xUpperBound_m2, xId2MemberId, memberForceByUnitWeight, memberLength_m, steelUnitWeight_knM3, steelYieldStrength_knM2 );

% LP + uncertainty in material property
steelYieldStrengthCov = 0.1;
[optimXArray_m2_uncertainMaterial, optimCostArray_uncertainMaterial, loadSampleArraySelected_uncertainMaterial, optimTime_uncertainMaterial, yieldStrengthSampleArray] = ...
    optimizeBpf_truss_linprog_uncertainMaterial( targetBpf, targetCov, activeAndFailRatio, solutionStepRatio, solutionChangeRatioTol, maxMemberForceInKNArray, xLowerBound_m2, xUpperBound_m2, xId2MemberId, memberForceByUnitWeight, memberLength_m, steelUnitWeight_knM3, steelYieldStrength_knM2, steelYieldStrengthCov );

% MIP + uncertainty in material property
xSolutionSet_m2 = 0.0010:0.0003:0.0040;
[optimXArray_m2_uncertainMaterial_milp, optimCostArray_uncertainMaterial_milp, loadSampleArraySelected_uncertainMaterial_milp, optimTime_uncertainMaterial_milp] = ...
    optimizeBpf_truss_milp_uncertainMaterial( targetBpf, targetCov, activeAndFailRatio, solutionChangeRatioTol, maxMemberForceInKNArray, xSolutionSet_m2, xId2MemberId, memberForceByUnitWeight, memberLength_m, steelUnitWeight_knM3, steelYieldStrength_knM2, steelYieldStrengthCov );


save ex3_optim