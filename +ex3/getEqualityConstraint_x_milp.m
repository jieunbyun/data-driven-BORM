function [Aeq_x, beq_x] = getEqualityConstraint_x_milp( nX, nSolution, nZ0, nZ, nS )

nDecisionVariable = nX*nSolution + nZ0 + nZ + nS;
Aeq_x = zeros( nX, nDecisionVariable );
beq_x = ones( nX, 1 );

for iXIndex = 1:nX
    iXLocationArray = nSolution * (iXIndex-1) + (1:nSolution );
    Aeq_x( iXIndex, iXLocationArray ) = 1;    
end
