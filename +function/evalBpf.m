function [bpf, z0, zArray, z0IndexInGValArray] = evalBpf( gValArray, nData )

if nargin < 2
    nData = length( gValArray );
end


[gValArraySort, gValArraySortIndex] = sort( gValArray, 'descend' );
gValArraySortCumSum = cumsum( gValArraySort );
gValArraySortCumSumFirstPositiveIndex = find( gValArraySortCumSum <= 0, 1 );

if gValArraySortCumSum( end ) > 0
    bpf = 1;
    z0 = -inf;
    zArray = inf( size(gValArray ) );
    z0IndexInGValArray = -inf;
    
elseif gValArraySortCumSum( 1 ) == 0
    bpf = 1/nData;
    z0 = 0;
    zArray = zeros( size( gValArray ) );
    z0IndexInGValArray = gValArraySortIndex( 1 );
    
elseif gValArraySortCumSum( 1 ) < 0
    bpf = 0;
    z0 = inf;
    zArray = -inf( size(gValArray ) );
    z0IndexInGValArray = 0;

else
    
    z0 = gValArraySort( gValArraySortCumSumFirstPositiveIndex );
    zArray = gValArray - z0; zArray( zArray < 0 ) = 0;
    bpf = sum( zArray ) / nData / (-z0);
    z0IndexInGValArray = gValArraySortIndex( gValArraySortCumSumFirstPositiveIndex );
    
end
