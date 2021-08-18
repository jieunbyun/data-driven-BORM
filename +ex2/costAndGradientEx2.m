function [cost, grad] = costAndGradientEx2( x, nLimitFun, targetBpf, penalty, nSample )

c1 = 6.74135e-5; c2 = 2.93585e-6; b2 = 3.556e2;

x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
cost = c1 * x1^2 * x2 + c2 * x3 * x4 * ( b2+x2 );

grad = zeros( 1,length(x) );
grad(1) = c1 * (2*x1) * x2;

grad(2) = c1 * x1^2 + c2 * x3 * x4;
grad(3) = c2 * x4 * ( b2+x2 );
grad(4) = c2 * x3 * ( b2+x2 );


if nargin > 3 % Penalized cost to find a feasible solution
    nX = 4;
    z0LocationArray = ( nX+(1:nLimitFun ) );
    zLocationArray = ( nX + nLimitFun + 1 ):length(x);
    
    z0 = x( z0LocationArray ); zArray = x( zLocationArray );
    cost = cost + penalty * sum( z0 + 1/targetBpf/nSample*sum( zArray ) );
    
    grad( z0LocationArray ) = penalty;
    grad( zLocationArray ) = penalty / targetBpf / nSample;
end

    