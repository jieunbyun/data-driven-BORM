function [funValArray, funGradArray] = limitFunEx2_fun3( vRowArray, gradientTrueOrFalse )

v1 = vRowArray(:,1); v4 = vRowArray(:,4);

funValArray = v1 ./ v4 - 1;

funGradArray = [];
if ( nargin > 1 ) && gradientTrueOrFalse
    nData = size( vRowArray, 1 );
    
    dG3V1 = 1./ v4; dG3V4 = - v1 ./ v4.^2; 
    funGradArray = [dG3V1 zeros(nData,1) zeros(nData,1) dG3V4];
    
end