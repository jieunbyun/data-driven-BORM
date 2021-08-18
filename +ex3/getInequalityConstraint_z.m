function [Aineq_z, bineq_z] = getInequalityConstraint_z( nX, nMember, nSampleActive, xId2MemberId, steelYieldStrength_knM2 )

nZ0 = nMember; nZ = nMember*nSampleActive; nS = nMember*nSampleActive;
nOptimVariable = nX + nZ0+nZ+nS;

Aineq_z = zeros( nZ, nOptimVariable ); bineq_z = zeros( nZ, 1 );
ijConstLocation = 0;
for iMemberIndex = 1:nMember
    iZ0Location = nX + iMemberIndex;
    iXIndex = cellfun( @(x) ismember( iMemberIndex,x ), xId2MemberId );
    
    for jSampleIndex = 1:nSampleActive
        ijConstLocation = ijConstLocation + 1;
        Aineq_z( ijConstLocation, iZ0Location ) = -1;
        Aineq_z( ijConstLocation, nX + nZ0 + nZ + ijConstLocation ) = 1;
        Aineq_z( ijConstLocation, nX + nZ0 + ijConstLocation ) = -1;
        Aineq_z( ijConstLocation, iXIndex ) = -steelYieldStrength_knM2;
    end
end