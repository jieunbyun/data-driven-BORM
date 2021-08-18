function [Aineq_s, bineq_s] = getInequalityConstraint_s( xId2MemberId, loadSampleArray, activeSampleArrayIndex, memberForceByUnitWeight, memberLengthArray, steelUnitWeight_knM3 )

nX = length( xId2MemberId );
nMember = size( loadSampleArray, 1 );
nSampleActive = size( activeSampleArrayIndex, 2 );


nZ0 = nMember; nZ = nMember*nSampleActive; nS = nMember*nSampleActive;
nOptimVariable = nX + nZ0+nZ+nS;

Aineq_s = zeros( 2*nS, nOptimVariable ); bineq_s = zeros( 2*nS, 1 );
ijConstLocation = 0;
ijSIndex = 0;
for iMemberIndex = 1:nMember
    
    iMemberForceByUnitWeight = memberForceByUnitWeight( iMemberIndex, : );
    iXWeight = zeros( 1,nX );
    for jXIndex = 1:nX
        jMemberArray = xId2MemberId{ jXIndex };
        ijWeight = steelUnitWeight_knM3 * sum( iMemberForceByUnitWeight( jMemberArray ) .* memberLengthArray( jMemberArray ) );
        iXWeight( jXIndex ) = ijWeight;
    end
    
    iActiveSampleIndex = activeSampleArrayIndex( iMemberIndex, : );
    for jSampleIndex = 1:nSampleActive
        ijSIndex = ijSIndex + 1;
        ijSLocation = nX + nZ0 + nZ + ijSIndex;
        
        ijConstLocation = ijConstLocation + 1;
        Aineq_s( ijConstLocation, 1:nX ) = iXWeight;
        
        Aineq_s( ijConstLocation, ijSLocation) = -1;
        
        ijActiveSampleIndex = iActiveSampleIndex( jSampleIndex );
        ijLoad = loadSampleArray( iMemberIndex, ijActiveSampleIndex );
        bineq_s( ijConstLocation ) = - ijLoad;
        
        ijConstLocation = ijConstLocation + 1;
        Aineq_s( ijConstLocation, 1:nX ) = -iXWeight;
        Aineq_s( ijConstLocation, ijSLocation) = -1;
        bineq_s( ijConstLocation ) = ijLoad;
        
    end
end
