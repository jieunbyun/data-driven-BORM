%{
13 May 21
Ji-Eun Byun

Ex2 - Sensitivity analysis
%}

clear;
import ex2.*
import function.*

load ex2_optim

dBpf_dSigma = zeros( nLimitFun, nV );
dBpf_dSigma_numerical = zeros( nLimitFun, nV );
dBpf_dX = zeros( nLimitFun, nX );
dBpf_dX_numerical = zeros( nLimitFun, nX );

%% Sensitivity analysis
vRowArray2 = repmat( vStdArray, size( USample2,1 ), 1 ) .* USample2 + repmat( optimXArray2(:).', size( USample2,1 ), 1 );
delta = 1e-6;

for kLimitFunIndex = 1:nLimitFun
    
    kLimitFunValArray2 = limitFunEx2_singleFun( vRowArray2, false, kLimitFunIndex );
    kBpf2 = evalBpf( kLimitFunValArray2 );
    
    if kBpf2

        % Sigma
        kDBpf_dSigma = zeros( 1,nV );
        for iSigmaIndex = 1:nV
            kiTheta = strcat( 'sigma', num2str(iSigmaIndex) );
            kiDBpf_dSigma = limitFunGradEx2( USample2, optimXArray2, vStdArray, kLimitFunIndex, kiTheta );
            kDBpf_dSigma( iSigmaIndex ) = kiDBpf_dSigma;
        end
        dBpf_dSigma( kLimitFunIndex, : ) = kDBpf_dSigma;


        kDBpf_dSigma_numerical = zeros( 1,nV );
        for iSigmaIndex = 1:nV
            kiDelta = zeros(1,nV); kiDelta( iSigmaIndex ) = delta;
            kiVStdArray_delta = vStdArray + kiDelta;
            kiVRowArray2_delta = repmat( kiVStdArray_delta, size( USample2,1 ), 1 ) .* USample2 + repmat( optimXArray2(:).', size( USample2,1 ), 1 );
            kiLimitFunValArray2_delta = limitFunEx2_singleFun(kiVRowArray2_delta,false,kLimitFunIndex);
            kiBpf2_delta = evalBpf( kiLimitFunValArray2_delta );
            kiDBpf_dSigma_numerical = (kiBpf2_delta - kBpf2) /delta;

            kDBpf_dSigma_numerical( iSigmaIndex ) = kiDBpf_dSigma_numerical;
        end
        dBpf_dSigma_numerical( kLimitFunIndex, : ) = kiDBpf_dSigma_numerical;


        disp( ['%%% Limit-state function ', num2str(kLimitFunIndex), ' w.r.t. Sigmas: Derivation, numerical, and relative difference %%%'] )
        disp( num2str( kDBpf_dSigma) )
        disp( num2str( kDBpf_dSigma_numerical) )
        disp( num2str( abs(kDBpf_dSigma-kDBpf_dSigma_numerical)./kDBpf_dSigma ) )


        % X
        kDBpf_dX_fun1 = zeros( 1,nV );
        for iSigmaIndex = 1:nV
            kiTheta = strcat( 'x', num2str(iSigmaIndex) );
            kiDBpf_dSigma = limitFunGradEx2( USample2, optimXArray2, vStdArray, kLimitFunIndex, kiTheta );
            kDBpf_dX_fun1( iSigmaIndex ) = kiDBpf_dSigma;
        end
        dBpf_dX( kLimitFunIndex,: ) = kiDBpf_dSigma;

        kDBpf_dX_fun1_numerical = zeros( 1,nV );
        for iSigmaIndex = 1:nV
            kiDelta = zeros(1,nV); kiDelta( iSigmaIndex ) = delta;
            iOptimXArray2 = optimXArray2(:).' + kiDelta;
            kiVRowArray2_delta = repmat( vStdArray, size( USample2,1 ), 1 ) .* USample2 + repmat( iOptimXArray2, size( USample2,1 ), 1 );
            kiLimitFunValArray2_delta = limitFunEx2_singleFun(kiVRowArray2_delta,false,kLimitFunIndex);
            kiBpf2_delta = evalBpf( kiLimitFunValArray2_delta );
            kiDBpf_dSigma_numerical = (kiBpf2_delta - kBpf2) /delta;

            kDBpf_dX_fun1_numerical( iSigmaIndex ) = kiDBpf_dSigma_numerical;
        end
        dBpf_dX_numerical( kLimitFunIndex,: ) = kDBpf_dX_fun1_numerical;

        disp( ['%%% Limit-state function ', num2str(kLimitFunIndex), ' w.r.t. Xs: Derivation, numerical, and relative difference %%%'] )
        disp( num2str( kDBpf_dX_fun1) )
        disp( num2str( kDBpf_dX_fun1_numerical) )
        disp( num2str( abs(kDBpf_dX_fun1-kDBpf_dX_fun1_numerical)./kDBpf_dX_fun1 ) )
    end
end

save ex2_sensitivity