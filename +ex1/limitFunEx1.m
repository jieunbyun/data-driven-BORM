function funValArray = limitFunEx1( vRowArray )

nV = 2;
if ~isempty( vRowArray ) && size( vRowArray, 2 ) ~= nV
    error( 'V Array must have two columns' )
end

import ex1.*
funVal1 = limitFunEx1_fun1( vRowArray, false);
funVal2 = limitFunEx1_fun2( vRowArray, false);
funValArray = [funVal1 funVal2];