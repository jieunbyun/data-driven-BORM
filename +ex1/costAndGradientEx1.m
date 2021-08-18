function [cost, grad] = costAndGradientEx1( x, nLimitFun, targetBpf, penalty, nSample )

x1 = x(1); x2 = x(2);

cost = (x1 - 3.7)^2 + (x2 - 4)^2;

grad = zeros(1, length(x) );
grad(1) = 2*(x1 - 3.7);
grad(2) = 2*(x2 - 4);


if nargin > 3 % Penalized cost to find a feasible solution
    nX = 2;
    z0LocationArray = ( nX+(1:nLimitFun ) );
    zLocationArray = ( nX + nLimitFun + 1 ):length(x);
        
    z0 = x(nX+1); zArray = x( (nX+1+1):end );
    cost = cost + penalty * ( z0 + 1/targetBpf/nSample*sum( zArray ) );
    
    grad( z0LocationArray ) = penalty;
    grad( zLocationArray ) = penalty / targetBpf / nSample;
end

    