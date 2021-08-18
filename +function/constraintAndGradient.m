function [constIneq, constEq, gradIneq, gradEq] = constraintAndGradient( x, uRowArray, activeSampleIndexArray, evalLimitFun_singleFun, nX, vStdArray, targetBpf, nSample )

nLimitFun = size( activeSampleIndexArray, 2 );   

constEq = [];
gradEq = [];

nDecisionVariable = length( x ); nActiveSample = size( activeSampleIndexArray, 1 );

nConstraint = nActiveSample*nLimitFun; 
constIneq = zeros( nConstraint, 1 );
gradIneq = zeros( nDecisionVariable, nConstraint );


for iLimitFunIndex = 1:nLimitFun
    
    iZ0 = x( nX + iLimitFunIndex );    
    iZLocationArray = nX + nLimitFun + nActiveSample * (iLimitFunIndex-1) + (1:nActiveSample );
    iZArray = x( iZLocationArray );
    
    iActiveSampleIndexArray = activeSampleIndexArray( :, iLimitFunIndex );
    for jZIndex = 1:nActiveSample
        ijActiveSampleIndex = iActiveSampleIndexArray( jZIndex );
        ijVArray = x(1:nX).' + vStdArray .* uRowArray( ijActiveSampleIndex,: );
        [ijFunVal, ijFunGrad] = evalLimitFun_singleFun( ijVArray, true, iLimitFunIndex );
        
        ijZ = iZArray( jZIndex );
        ijConst = ijFunVal - iZ0 - ijZ;
        
        ijConstLocation = nActiveSample*(iLimitFunIndex-1) + jZIndex;
        constIneq( ijConstLocation ) = ijConst;
        
        gradIneq( 1:nX, ijConstLocation ) = ijFunGrad;
        gradIneq( nX+iLimitFunIndex, ijConstLocation ) = -1;
        ijZLocation = iZLocationArray( jZIndex );
        gradIneq( ijZLocation, ijConstLocation ) = -1;
    end
end

if nargin > 6 % Relaxed reliability constraint
    
    constIneq_penalty = zeros( nLimitFun, 1 );
    gradIneq_penalty = zeros( nDecisionVariable, nLimitFun );
    
    for iLimitFunIndex = 1:nLimitFun
    
        iZ0 = x( nX + iLimitFunIndex );    
        iZLocationArray = nX + nLimitFun + nActiveSample * (iLimitFunIndex-1) + (1:nActiveSample );
        iZArray = x( iZLocationArray );
        iY = x( nDecisionVariable-(nLimitFun-iLimitFunIndex) );
        
        iConst = iZ0 + 1/targetBpf/nSample*sum(iZArray) - iY;
        
        
        iGradIneq_penalty = zeros( nDecisionVariable, 1 );
        iGradIneq_penalty( nX + iLimitFunIndex ) = 1;
        iGradIneq_penalty( iZLocationArray ) = 1/targetBpf/nSample;
        iGradIneq_penalty( nDecisionVariable-(nLimitFun-iLimitFunIndex) ) = -1;
        
        constIneq_penalty( iLimitFunIndex ) = iConst;
        gradIneq_penalty( :,iLimitFunIndex ) = iGradIneq_penalty;
        
    end
    
    constIneq = [constIneq_penalty; constIneq];
    gradIneq = [gradIneq_penalty gradIneq];
end


