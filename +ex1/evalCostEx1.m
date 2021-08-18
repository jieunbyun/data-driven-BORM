function cost = evalCostEx1( xArray )

nX = 2;
if length( xArray ) ~= nX
    error( 'Given x vector must have 2 elements, i.e. the number of decision variables' )
end

cost = (xArray(1) - 3.7)^2 + (xArray(2) - 4)^2;