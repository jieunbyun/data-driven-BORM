function [funValArray, funGradArray] = limitFunEx2_fun2( vRowArray, gradientTrueOrFalse )

v3 = vRowArray(:,3); v4 = vRowArray(:,4);

b1 = 2.6688e4; b2 = 3.556e2; b7 = 2.0685e2;

sigma = 6 * b1 * b2 ./ v3.^2 ./v4;
funValArray = sigma / b7 - 1;


funGradArray = [];
if ( nargin > 1 ) && gradientTrueOrFalse
    nData = size( vRowArray, 1 );
    
    dSigmaV3 = -2 * 6 * b1 * b2 ./ v3.^3 ./v4;
    dSigmaV4 = - 6 * b1 * b2 ./ v3.^2 ./ v4.^2;
    dSigma = [zeros(nData,1) zeros(nData,1) dSigmaV3 dSigmaV4];

    funGradArray = dSigma / b7;
end