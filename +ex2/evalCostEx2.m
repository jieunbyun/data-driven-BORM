function cost = evalCostEx2( xArray )

nX = 4;
if length( xArray ) ~= nX
    error( 'Given x vector must have 2 elements, i.e. the number of decision variables' )
end

c1 = 6.74135e-5; c2 = 2.93585e-6; b2 = 3.556e2;
x1 = xArray(1); x2 = xArray(2); x3 = xArray(3); x4 = xArray(4);
cost = c1 * x1^2 * x2 + c2 * x3 * x4 * ( b2+x2 );