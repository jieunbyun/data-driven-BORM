function [funValArray, funGradArray] = limitFunEx1_fun2( vRowArray, gradientTrueOrFalse )

nV = 2;
v1 = vRowArray(:,1); v2 = vRowArray(:,2);
funValArray = - v1 - v2 + 3;

funGradArray = [];
if gradientTrueOrFalse
    nSample = size( vRowArray, 1 );
    funGradArray = -ones( nSample, nV );
end
