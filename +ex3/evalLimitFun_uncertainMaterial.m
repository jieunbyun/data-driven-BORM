function limitStateFunValArray = evalLimitFun_uncertainMaterial( xValArray, xId2MemberId, loadSampleArray, memberForceByUnitWeight, memberLengthArray, steelUnitWeight_knM3, yieldStrengthSampleArray )

nX = length( xValArray );
nMember = size( loadSampleArray, 1 );
nSample = size( loadSampleArray, 2 );

memberAreaArray = zeros( 1,nMember);
for iXIndex = 1:nX
    iMemberArray = xId2MemberId{ iXIndex };
    iX = xValArray( iXIndex );
    memberAreaArray( iMemberArray ) = iX;
end

limitStateFunValArray = zeros( nMember, nSample );
for iSampleIndex = 1:nSample
    iLoad = loadSampleArray( :, iSampleIndex );
    for jMemberIndex = 1:nMember
        ijDemandByWeight = steelUnitWeight_knM3 * sum( memberForceByUnitWeight( jMemberIndex, : ) .* memberAreaArray .* memberLengthArray );
        ijDemandByTraffic = iLoad( jMemberIndex );
        
        ijYieldStrength = yieldStrengthSampleArray( jMemberIndex, iSampleIndex );
        ijLimitFunVal = abs( ijDemandByWeight + ijDemandByTraffic ) - ijYieldStrength * memberAreaArray( jMemberIndex );
        limitStateFunValArray( jMemberIndex, iSampleIndex ) = ijLimitFunVal;
    end
end