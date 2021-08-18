function [limitStateFunValArray, demandArray] = evalLimitFun( xValArray, xId2MemberId, loadSampleArray, memberForceByUnitWeight, memberLengthArray, steelUnitWeight_knM3, steelYieldStrength_knM2 )

nX = length( xValArray );
nMember = size( loadSampleArray, 1 );
nSample = size( loadSampleArray, 2 );

if isscalar(steelYieldStrength_knM2)
    steelYieldStrength_knM2 = repmat( steelYieldStrength_knM2, nMember, nSample );
end


memberAreaArray = zeros( 1,nMember);
for iXIndex = 1:nX
    iMemberArray = xId2MemberId{ iXIndex };
    iX = xValArray( iXIndex );
    memberAreaArray( iMemberArray ) = iX;
end

limitStateFunValArray = zeros( nMember, nSample ); demandArray = zeros( nMember, nSample );
for iSampleIndex = 1:nSample
    iLoad = loadSampleArray( :, iSampleIndex );
    for jMemberIndex = 1:nMember
        ijDemandByWeight = steelUnitWeight_knM3 * sum( memberForceByUnitWeight( jMemberIndex, : ) .* memberAreaArray .* memberLengthArray );
        ijDemandByTraffic = iLoad( jMemberIndex );
        ijSteelYieldStrength = steelYieldStrength_knM2( jMemberIndex, iSampleIndex );
        
        ijDemand = ijDemandByWeight + ijDemandByTraffic;
        ijLimitFunVal = abs( ijDemand ) - ijSteelYieldStrength * memberAreaArray( jMemberIndex );
        
        demandArray( jMemberIndex, iSampleIndex ) = ijDemand;
        limitStateFunValArray( jMemberIndex, iSampleIndex ) = ijLimitFunVal;
    end
end