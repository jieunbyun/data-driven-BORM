%{
21 Apr 21
Ji-Eun Byun

Ex1 - Sensitivity analysis
%}

clear;
import ex1.*
import function.*

load ex1_optim

%% Sensitivity analysis
limitFunIndex = 1;
vRowArray2 = repmat( vStdArray, size( USample2,limitFunIndex ), 1 ) .* USample2 + repmat( optimXArray2(:).', size( USample2,limitFunIndex ), 1 );
limitFunValArray2 = limitFunEx1_fun1( vRowArray2 );
bpf2 = evalBpf( limitFunValArray2 );

% Sigma
dBpf_dSigma1_fun1 = limitFunGradEx1_fun1( USample2, optimXArray2, vStdArray, 'sigma1' );
delta = 1e-6;
vStdArray_delta = vStdArray + delta * [1 0];
vRowArray2_delta = repmat( vStdArray_delta, size( USample2,limitFunIndex ), 1 ) .* USample2 + repmat( optimXArray2(:).', size( USample2,limitFunIndex ), 1 );
limitFunValArray2_delta = limitFunEx1_fun1( vRowArray2_delta );
bpf2_delta = evalBpf( limitFunValArray2_delta );
dBpf_dSigma1_fun1_numerical = (bpf2_delta - bpf2) /delta;

dBpf_dSigma2_fun1 = limitFunGradEx1_fun1( USample2, optimXArray2, vStdArray, 'sigma2' );
vStdArray_delta = vStdArray + delta * [0 1];
vRowArray2_delta = repmat( vStdArray_delta, size( USample2,limitFunIndex ), 1 ) .* USample2 + repmat( optimXArray2(:).', size( USample2,limitFunIndex ), 1 );
limitFunValArray2_delta = limitFunEx1_fun1( vRowArray2_delta );
bpf2_delta = evalBpf( limitFunValArray2_delta );
dBpf_dSigma2_fun1_numerical = (bpf2_delta - bpf2) /delta;

disp( '%%% Limit-state function 1: Derivation, numerical, and relative difference %%%' )
disp( ['Sigma 1: ' num2str( [dBpf_dSigma1_fun1 dBpf_dSigma1_fun1_numerical abs(dBpf_dSigma1_fun1-dBpf_dSigma1_fun1_numerical)/dBpf_dSigma1_fun1_numerical] ) ] )
disp( ['Sigma 2: ' num2str( [dBpf_dSigma2_fun1 dBpf_dSigma2_fun1_numerical abs(dBpf_dSigma2_fun1-dBpf_dSigma2_fun1_numerical)/dBpf_dSigma2_fun1_numerical] ) ]  )

% X
dBpf_dX1_fun1 = limitFunGradEx1_fun1( USample2, optimXArray2, vStdArray, 'x1' );
xArray_delta = optimXArray2 + delta * [1 0].';
vRowArray2_delta = repmat( vStdArray, size( USample2,limitFunIndex ), 1 ) .* USample2 + repmat( xArray_delta(:).', size( USample2,limitFunIndex ), 1 );
limitFunValArray2_delta = limitFunEx1_fun1( vRowArray2_delta );
bpf2_delta = evalBpf( limitFunValArray2_delta );
dBpf_dX1_fun1_numerical = (bpf2_delta - bpf2) /delta;

dBpf_dX2_fun1 = limitFunGradEx1_fun1( USample2, optimXArray2, vStdArray, 'x2' );
XArray_delta = optimXArray2 + delta * [0 1].';
vRowArray2_delta = repmat( vStdArray, size( USample2,limitFunIndex ), 1 ) .* USample2 + repmat( XArray_delta(:).', size( USample2,limitFunIndex ), 1 );
limitFunValArray2_delta = limitFunEx1_fun1( vRowArray2_delta );
bpf2_delta = evalBpf( limitFunValArray2_delta );
dBpf_dX2_fun1_numerical = (bpf2_delta - bpf2) /delta;

disp( '%%% Limit-state function 1: Derivation and Numerical %%%' )
disp( ['X 1: ' num2str( [dBpf_dX1_fun1 dBpf_dX1_fun1_numerical abs(dBpf_dX1_fun1-dBpf_dX1_fun1_numerical)/dBpf_dX1_fun1] ) ] )
disp( ['X 2: ' num2str( [dBpf_dX2_fun1 dBpf_dX2_fun1_numerical abs(dBpf_dX2_fun1-dBpf_dX2_fun1_numerical)/dBpf_dX2_fun1] ) ]  )

save ex1_sensitivity