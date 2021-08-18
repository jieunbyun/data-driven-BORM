function [optimPf, optimBpf, optimLimitFunValArray, z0, zArray, z0Index, vArray] = summarizeOptimResult( optimXArray, USample, optimTime, funEvalCost, evalLimitFun, vStdArray, targetPf, targetBpf, optimRoundName )

import function.*

optimX = optimXArray( :, end).';             
optimX = round( optimX, 2 );
nSample = size( USample, 1 );                                                                              
vArray = repmat( vStdArray, nSample, 1 ) .* USample + repmat( optimX, nSample, 1 );
optimLimitFunValArray = evalLimitFun( vArray );
optimPf = mean( optimLimitFunValArray > 0 );

nLimitFun = size( optimLimitFunValArray, 2 );
optimBpf = zeros( 1, nLimitFun ); z0 = zeros( 1,nLimitFun ); zArray = zeros( nSample, nLimitFun ); z0Index = zeros( 1,nLimitFun );  
for iLimitFunIndex = 1:nLimitFun
    iFunValArray = optimLimitFunValArray( :, iLimitFunIndex );
    [iOptimBpf, iZ0, iZArray, iZ0Index] = evalBpf( iFunValArray );
    optimBpf( iLimitFunIndex ) = iOptimBpf;
    z0( iLimitFunIndex ) = iZ0;
    zArray( :, iLimitFunIndex ) = iZArray;
    z0Index( iLimitFunIndex ) = iZ0Index;
end
cost = funEvalCost( optimX );

disp( ['%%%%%% Optimization: ' optimRoundName ' %%%%%%' ] )
disp( ['Number of iterations: ' num2str( size( optimXArray,2 ) ) ' / Optimization time (sec): ' num2str( optimTime )] )
disp( ['Optimal X: ' num2str( optimX ) ' / Cost: ' num2str( cost )] )
disp( ['Pf: ' num2str( optimPf ) ' (target Pf: ' num2str( targetPf ) ') / bPf: ' num2str( optimBpf ) ' (target bPf: ' num2str( targetBpf ) ') / z0: ' num2str( z0 )] )
disp( ['tau: ' num2str( optimBpf ./ optimPf )] )