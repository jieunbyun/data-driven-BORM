function tau = getBufferedTailIndexFromPf( pf )

if pf < 1e-6 || pf > 0.5
    error( 'Given pf must lie in [1e-6 5e-1]' )
else
    tau = interp1( log([1e-6 0.01 0.3 0.5]), [2.68 2.61 2.4 2], log(pf) );
end
    