function memberForceByUnitWeight = evalMemberForceByUnitWeight( memberConnectivity, memberForceArray )

nMember = size( memberConnectivity, 1 );
memberForceByUnitWeight = zeros( nMember, nMember );
for iWeightMemberIndex = 1:nMember
    iMemberEndNodePair = memberConnectivity( iWeightMemberIndex, : );
    
    for jAffectedMemberIndex = 1:nMember
        ijMemberForce = sum( memberForceArray( iMemberEndNodePair, jAffectedMemberIndex ) );
        memberForceByUnitWeight( jAffectedMemberIndex, iWeightMemberIndex ) = ijMemberForce;
    end
end


