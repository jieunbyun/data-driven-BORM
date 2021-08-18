function costFun = getCostFunction_milp( nSampleActive, xId2MemberId, memberLengthArray, steelUnitWeight_knM3, xSolutionSet )

nX = length( xId2MemberId );
nMember = length( memberLengthArray );
nArea = length( xSolutionSet );

costFun_x = zeros( 1, nX*nArea );
iXLocation = 0;
for iXIndex = 1:nX
    for jAreaIndex = 1:nArea
        iMemberArray = xId2MemberId{ iXIndex };
        iWeight = steelUnitWeight_knM3 * sum( memberLengthArray( iMemberArray ) );
        jArea = xSolutionSet( jAreaIndex );
        
        iXLocation = iXLocation + 1;
        costFun_x( iXLocation ) = iWeight*jArea;
    end
end

nZ0 = nMember; nZ = nMember*nSampleActive; nS = nMember*nSampleActive;
costFun = [costFun_x zeros( 1,nZ0+nZ+nS )];