function funValArray = limitFunEx2( vRowArray )

nV = 4;
if ~isempty( vRowArray ) && size( vRowArray, 2 ) ~= nV
    error( 'V Array must have two columns' )
end

import ex2.*
funVal1 = limitFunEx2_fun1( vRowArray, false );
funVal2 = limitFunEx2_fun2( vRowArray, false );
funVal3 = limitFunEx2_fun3( vRowArray, false );
funVal4 = limitFunEx2_fun4( vRowArray, false );
funVal5 = limitFunEx2_fun5( vRowArray, false );
funValArray = [funVal1 funVal2 funVal3 funVal4 funVal5];