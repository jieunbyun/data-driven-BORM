function z0AndZ = evalZValues( gValArray, nData, bPf )

gValArray = gValArray(:);

nGVal = length( gValArray );
gValArrayFull = [gValArray; ( min( gValArray ) - 1 ) * ones( nData-nGVal, 1 )];

z0 = quantile( gValArrayFull, 1-bPf );
zArray = gValArray - z0;
zArray( zArray<0 ) = 0;

z0AndZ = [z0; zArray];