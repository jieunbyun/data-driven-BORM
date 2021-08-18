function dBpf_dTheta = limitFunGradEx2( uRowArray, x, vStdArray, limitFunIndex, theta )

import ex2.*
import function.*

nData = size( uRowArray, 1 );
vStdArray = vStdArray(:).'; x = x(:).';
vRowArray = repmat( vStdArray, nData, 1 ) .* uRowArray + repmat( x, nData, 1 );

[funValArray, funGradArray] = limitFunEx2_singleFun( vRowArray, true, limitFunIndex );
[~, z0, zArray, z0Index] = evalBpf( funValArray );


switch theta
    case 'sigma1'
        dG_dThetaArray = funGradArray(:,1) .* uRowArray(:,1);
    case 'sigma2'
        dG_dThetaArray = funGradArray(:,2) .* uRowArray(:,2);
    case 'sigma3'
        dG_dThetaArray = funGradArray(:,3) .* uRowArray(:,3);    
    case 'sigma4'
        dG_dThetaArray = funGradArray(:,4) .* uRowArray(:,4);
    
    case 'x1'
        dG_dThetaArray = funGradArray(:,1);
    case 'x2'
        dG_dThetaArray = funGradArray(:,2);
    case 'x3'
        dG_dThetaArray = funGradArray(:,3);
    case 'x4'
        dG_dThetaArray = funGradArray(:,4);
        
    otherwise
        error( 'Given parameter is not in database' )
end
        
zPositiveIndex = ( zArray > 0 );
dBpf_dTheta = 1/nData * sum( dG_dThetaArray( zPositiveIndex ) / (-z0) ) + 1/nData * dG_dThetaArray( z0Index ) * sum( funValArray( zPositiveIndex ) ) / z0^2;