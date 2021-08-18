function [funValArray, funGradArray] = limitFunEx1_singleFun( vRowArray, gradientTrueOrFalse, limitFunIndex )

nV = 2;
if ~isempty( vRowArray ) && size( vRowArray, 2 ) ~= nV
    error( 'V Array must have two columns' )
end

nLimitFun = 2;
if ~isscalar( limitFunIndex )
    error( 'Limit-state function index must be scalar' )
elseif ~ismember( limitFunIndex, 1:nLimitFun )
    error( 'Given limit-state function index is not in database' )
end

if nargin < 3
    gradientTrueOrFalse = false;
end

import ex1.*
switch limitFunIndex
    case 1
        [funValArray, funGradArray] = limitFunEx1_fun1( vRowArray, gradientTrueOrFalse );
    case 2
        [funValArray, funGradArray] = limitFunEx1_fun2( vRowArray, gradientTrueOrFalse );
end