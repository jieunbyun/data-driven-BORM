function dBpf_dTheta = limitFunGradEx3( loadSampleArray, xValArray, memberIndex, theta, xId2MemberId, memberForceByUnitWeight, memberLengthArray, steelUnitWeight_knM3, steelYieldStrength_knM2, yieldStrengthSampleArray, steelYieldStrengthCov, steelYieldStrengthMean, demandLoadArray )

import ex3.*
import function.*

nData = size( loadSampleArray, 2 );

funValArray = evalLimitFun( xValArray, xId2MemberId, loadSampleArray, memberForceByUnitWeight, memberLengthArray, steelUnitWeight_knM3, yieldStrengthSampleArray );
funValArray = funValArray( memberIndex,: );
[~, z0, zArray, z0Index] = evalBpf( funValArray(:) );

xID = find( cellfun( @(x) ismember(memberIndex,x), xId2MemberId ) );

if strcmp( theta(1), 'x' )
    paramXID = str2num( theta(2:end) );
    theta = 'x';
end

switch theta
    case 'sigma'
        uSampleArray = ( yieldStrengthSampleArray( memberIndex,: ) - steelYieldStrength_knM2 ) / steelYieldStrength_knM2 / steelYieldStrengthCov;
        uSampleArray = uSampleArray(:);
        
        dG_dThetaArray = -xValArray( xID ) * steelYieldStrengthMean * uSampleArray;
        
    case 'x'
        memberIDArray = xId2MemberId{ paramXID };
        dG_dx = steelUnitWeight_knM3 * memberLengthArray( memberIDArray ) .* memberForceByUnitWeight( memberIndex,memberIDArray );
        dG_dx = sum( dG_dx ) * sign( demandLoadArray( memberIndex,: ) );
        
        if ismember( memberIndex, memberIDArray )
            dG_dThetaArray = dG_dx - yieldStrengthSampleArray( memberIndex, : );
        else
            dG_dThetaArray = dG_dx;
        end
        dG_dThetaArray = dG_dThetaArray(:);
        
    otherwise
        error( 'Given parameter is not in database' )
end
        
zPositiveIndex = ( zArray > 0 );
dBpf_dTheta = 1/nData * sum( dG_dThetaArray( zPositiveIndex ) / (-z0) ) + 1/nData * dG_dThetaArray( z0Index ) * sum( funValArray( zPositiveIndex ) ) / z0^2;