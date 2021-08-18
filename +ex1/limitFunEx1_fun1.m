function [funValArray, funGradArray] = limitFunEx1_fun1( vRowArray, gradientTrueOrFalse )

v1 = vRowArray(:,1); v2 = vRowArray(:,2);
funValArray = v1 .* sin( 4*v1 ) + 1.1*v2.*sin( 2*v2 );

funGradArray = [];
if nargin > 1
    if gradientTrueOrFalse
        dFun_V1 = sin( 4*v1 ) + v1 .* cos( 4*v1 ) *4;
        dFun_V2 = 1.1*sin( 2*v2 ) + 1.1*v2.*cos( 2*v2 )*2;

        funGradArray = [dFun_V1 dFun_V2];
    end
end
