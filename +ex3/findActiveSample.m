function [activeSampleArrayVal, activeSampleArrayIndex] = findActiveSample( limitStateFunValArray, nSampleActive )

nMember = size( limitStateFunValArray, 1 );

activeSampleArrayVal = zeros( nMember, nSampleActive );
activeSampleArrayIndex = zeros( nMember, nSampleActive );

for iMemberIndex = 1:nMember
    iLimitStateFunValArray = limitStateFunValArray( iMemberIndex, : );
    [iLimitStateFunValArraySort, iLimitStateFunValArraySortIndex] = sort( iLimitStateFunValArray, 'descend');
    
    activeSampleArrayVal( iMemberIndex, : ) = iLimitStateFunValArraySort( 1:nSampleActive );
    activeSampleArrayIndex( iMemberIndex, : ) = iLimitStateFunValArraySortIndex( 1:nSampleActive );
end