function memberLength = evalMemberLength( memberConnectivity, nodeCoord )

nMember = size( memberConnectivity, 1 );
memberLength = zeros( 1, nMember );
for iMemberIndex = 1:nMember
    iMemberEndNodePair = memberConnectivity( iMemberIndex, : );
    memberLength( iMemberIndex ) = norm( nodeCoord( iMemberEndNodePair(1), : ) - nodeCoord( iMemberEndNodePair(2),: ) );
end