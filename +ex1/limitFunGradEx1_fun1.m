function dBpf_dTheta = limitFunGradEx1_fun1( uRowArray, x, vStdArray, theta )

import ex1.*
import function.*

nData = size( uRowArray, 1 );
vStdArray = vStdArray(:).'; x = x(:).';
vRowArray = repmat( vStdArray, nData, 1 ) .* uRowArray + repmat( x, nData, 1 );

[funValArray, funGradArray] = limitFunEx1_fun1( vRowArray, true );
[~, z0, zArray, z0Index] = evalBpf( funValArray );

switch theta
    case 'sigma1'
        dG_dThetaArray = funGradArray(:,1) .* uRowArray(:,1);

    case 'sigma2'
        dG_dThetaArray = funGradArray(:,2) .* uRowArray(:,2);
        
    case 'x1'
        dG_dThetaArray = funGradArray(:,1);
        
    case 'x2'
        dG_dThetaArray = funGradArray(:,2);
        
    otherwise 
        error( 'Given parameter is not in database' )
        
end
zPositiveIndex = ( zArray > 0 );
dBpf_dTheta = 1/nData * sum( dG_dThetaArray( zPositiveIndex ) / (-z0) ) + 1/nData * dG_dThetaArray( z0Index ) * sum( funValArray( zPositiveIndex ) ) / z0^2;
