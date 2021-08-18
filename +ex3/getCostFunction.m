function costFun = getCostFunction( nSampleActive, xId2MemberId, memberLengthArray, steelUnitWeight_knM3 )

nX = length( xId2MemberId );
nMember = length( memberLengthArray );

costFun_x = zeros( 1, nX );
for iXIndex = 1:nX
    iMemberArray = xId2MemberId{ iXIndex };
    iWeight = steelUnitWeight_knM3 * sum( memberLengthArray( iMemberArray ) );
    costFun_x( iXIndex ) = iWeight;
end

nZ0 = nMember; nZ = nMember*nSampleActive; nS = nMember*nSampleActive;
costFun = [costFun_x zeros( 1,nZ0+nZ+nS )];