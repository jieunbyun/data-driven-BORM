function [funValArray, funGradArray] = limitFunEx2_singleFun( vRowArray, gradientTrueOrFalse, limitFunIndex )

nV = 4;
if ~isempty( vRowArray ) && size( vRowArray, 2 ) ~= nV
    error( 'V Array must have two columns' )
end

if ~isscalar( limitFunIndex )
    error( 'Limit-state function index must be scalar' )
end

if nargin < 3
    gradientTrueOrFalse = false;
end


import ex2.*
switch limitFunIndex
    case 1
        [funValArray, funGradArray] = limitFunEx2_fun1( vRowArray, gradientTrueOrFalse );
    case 2
        [funValArray, funGradArray] = limitFunEx2_fun2( vRowArray, gradientTrueOrFalse );
    case 3
        [funValArray, funGradArray] = limitFunEx2_fun3( vRowArray, gradientTrueOrFalse );
    case 4
        [funValArray, funGradArray] = limitFunEx2_fun4( vRowArray, gradientTrueOrFalse );
    case 5
        [funValArray, funGradArray] = limitFunEx2_fun5( vRowArray, gradientTrueOrFalse );
    otherwise
        error( 'Given limit-state function is not in database' )
end