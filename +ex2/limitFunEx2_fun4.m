function [funValArray, funGradArray] = limitFunEx2_fun4( vRowArray, gradientTrueOrFalse )

v3 = vRowArray(:,3); v4 = vRowArray(:,4);

b1 = 2.6688e4; b2 = 3.556e2; b3 = 2.0685e5;  b5 = 6.35; 

delta = 4 * b1 * b2^3 / b3 ./ v3.^3 ./ v4;

funValArray = delta / b5 - 1;


funGradArray = [];
if ( nargin > 1 ) && gradientTrueOrFalse
    nData = size( vRowArray, 1 );
    
    dDeltaV3 = -3 * 4 * b1 * b2^3 / b3 ./ v3.^4 ./ v4;
    dDeltaV4 = - 4 * b1 * b2^3 / b3 ./ v3.^3 ./ v4.^2;
    dDelta = [zeros(nData,1) zeros(nData,1) dDeltaV3 dDeltaV4];
        
    funGradArray = dDelta / b5;
    
end