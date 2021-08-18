function [Aineq_z0, bineq_z0] = getInequalityConstraint_z0( nDecisionVariable, nMember, nSample, nSampleActive, targetBpf )

nZ0 = nMember; nZ = nMember*nSampleActive; nS = nMember*nSampleActive;
nOptimVariable = nDecisionVariable + nZ0+nZ+nS;

Aineq_z0 = zeros( nZ0, nOptimVariable ); bineq_z0 = zeros( nZ0, 1 );
for iZ0Index = 1:nZ0
    iZ0Location = nDecisionVariable + iZ0Index;
    iZLocation = nDecisionVariable + nZ0 + nSampleActive * (iZ0Index-1 ) + (1:nSampleActive);
    
    Aineq_z0( iZ0Index, iZ0Location ) = 1;
    Aineq_z0( iZ0Index, iZLocation ) = 1/targetBpf/nSample;
end