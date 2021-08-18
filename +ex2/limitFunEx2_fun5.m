function [funValArray, funGradArray] = limitFunEx2_fun5( vRowArray, gradientTrueOrFalse )

v3 = vRowArray(:,3); v4 = vRowArray(:,4);

b1 = 2.6688e4; b2 = 3.556e2; b3 = 2.0685e5; b4 = 8.274e4;

Pc = 4.013 * v3 .* v4.^3 * sqrt(b3*b4) / 6 /b2^2 .* ( 1 - v3 / 4 / b2 * sqrt( b3/b4 ) );
funValArray = 1 - Pc / b1;


funGradArray = [];
if ( nargin > 1 ) && gradientTrueOrFalse
    nData = size( vRowArray, 1 );
    
    dPcV3 = 4.013 * v4.^3 * sqrt(b3*b4) / 6 /b2^2 .* ( 1 - v3 / 4 / b2 * sqrt( b3/b4 ) ) + ...
                   4.013 * v3 .* v4.^3 * sqrt(b3*b4) / 6 /b2^2 .* ( - 1 / 4 / b2 * sqrt( b3/b4 ) );
    dPcV4 = 4.013 * v3 .* (3*v4.^2) * sqrt(b3*b4) / 6 /b2^2 .* ( 1 - v3 / 4 / b2 * sqrt( b3/b4 ) );
    dPc = [zeros(nData,1) zeros(nData,1) dPcV3 dPcV4];
    
    funGradArray = - dPc / b1;
   
end
    
    
    