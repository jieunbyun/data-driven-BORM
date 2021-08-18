function [Aineq_z, bineq_z] = getInequalityConstraint_z_uncertainMaterial_milp(  xId2MemberId, yieldStrengthSampleArray, activeSampleArrayIndex, xSolutionSet )

nX = length( xId2MemberId )*length( xSolutionSet );
nMember = size( activeSampleArrayIndex, 1 );
nSampleActive = size( activeSampleArrayIndex, 2 );
nSolution = length( xSolutionSet );

nZ0 = nMember; nZ = nMember*nSampleActive; nS = nMember*nSampleActive;
nOptimVariable = nX + nZ0+nZ+nS;

Aineq_z = zeros( nZ, nOptimVariable ); bineq_z = zeros( nZ, 1 );
ijConstLocation = 0;
for iMemberIndex = 1:nMember
    iZ0Location = nX + iMemberIndex;
    iXIndex = find( cellfun( @(x) ismember( iMemberIndex,x ), xId2MemberId ) );
    
    for jSampleIndex = 1:nSampleActive
        ijConstLocation = ijConstLocation + 1;
        Aineq_z( ijConstLocation, iZ0Location ) = -1;
        Aineq_z( ijConstLocation, nX + nZ0 + nZ + ijConstLocation ) = 1;
        Aineq_z( ijConstLocation, nX + nZ0 + ijConstLocation ) = -1;
        
        ijActiveSampleArray = activeSampleArrayIndex( iMemberIndex, jSampleIndex );
        
        for kSolutionIndex = 1:nSolution
            ikXLocation = nSolution * (iXIndex-1) + kSolutionIndex;
            kSolution = xSolutionSet( kSolutionIndex );
            Aineq_z( ijConstLocation, ikXLocation ) = - kSolution * yieldStrengthSampleArray( iMemberIndex, ijActiveSampleArray );
        end
    end
end