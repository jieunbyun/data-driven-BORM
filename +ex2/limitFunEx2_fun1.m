function [funValArray, funGradArray] = limitFunEx2_fun1( vRowArray, gradientTrueOrFalse )

nData = size(vRowArray,1);
v1 = vRowArray(:,1); v2 = vRowArray(:,2); v3 = vRowArray(:,3);
b1 = 2.6688e4; b2 = 3.556e2; b6 = 9.377e1;
        
J = sqrt(2) * v1 .* v2 .* ( v2.^2/12 + ( v1 + v3 ).^2/4 );
M = b1 * ( b2 + v2/2 );
R = sqrt( v2.^2 + ( v1 + v3 ).^2 ) / 2;
t = b1 / sqrt(2) ./ v1 ./ v2;
tt = M .* R ./ J;
tau = sqrt( t.^2 + 2 * t .* tt .*v2 / 2 ./ R + tt.^2 );

funValArray =  tau / b6 - 1;

funGradArray = [];
if ( nargin > 1 ) && gradientTrueOrFalse

    dJV1 = sqrt(2) * v2 .* ( v2.^2/12 + ( v1 + v3 ).^2/4 ) + sqrt(2) * v1 .* v2 .* ( 2*( v1 + v3 )/4 );
    dJV2 = sqrt(2) * v1 .* ( v2.^2/12 + ( v1 + v3 ).^2/4 ) + sqrt(2) * v1 .* v2 .* ( (2*v2)/12 );
    dJV3 = sqrt(2) * v1 .* v2 .* ( 2*( v1 + v3 )/4 );
    dJ = [dJV1 dJV2 dJV3 zeros(nData,1)];

    dMV2 = b1 * ( 1/2 );
    dM = [0 dMV2 0 0];

    dRV1 = 1/2 * ( v2.^2 + ( v1 + v3 ).^2 ).^(-1/2) / 2 .* ( 2*(v1 + v3) );
    dRV2 = 1/2 * ( v2.^2 + ( v1 + v3 ).^2 ).^(-1/2) / 2 .* ( 2*v2 );
    dRV3 = dRV1;
    dR = [dRV1 dRV2 dRV3 zeros(nData,1)];

    dtV1 = - b1 / sqrt(2) ./ (v1.^2) ./ v2;
    dtV2 = - b1 / sqrt(2) ./ v1 ./ (v2.^2);
    dt = [dtV1 dtV2 zeros(nData,1) zeros(nData,1)];

    nData = size( vRowArray, 1 ); nV = 4;
    dtt = zeros( nData, nV );
    for iVIndex = 1:nV
        dttVi = dM( :, iVIndex ) .* R ./ J + dR( :, iVIndex ) .* M ./ J - dJ( :, iVIndex ) .* M .* R ./ J.^2;
        dtt( :, iVIndex ) = dttVi;
    end

    dtau = zeros( nData, nV );
    for iVIndex = 1:nV
        dtVi = dt(:,iVIndex); dttVi = dtt(:,iVIndex); dRVi = dR(:,iVIndex);
        if iVIndex ~= 2
            dtauVi =  (1/2) * ( t.^2 + 2 * t .* tt .*v2 / 2 ./ R + tt.^2 ).^(-1/2) .* ...
                            ( 2*t.*dtVi + 2 * dtVi .* tt .*v2 / 2 ./ R + 2 * t .* dttVi .*v2 / 2 ./ R - 2 * t .* tt .*v2 / 2 ./ R.^2 .* dRVi + (2*tt).*dttVi );
        else
            dtauVi =  (1/2) * ( t.^2 + 2 * t .* tt .*v2 / 2 ./ R + tt.^2 ).^(-1/2) .* ...
                            ( 2*t.*dtVi + 2 * dtVi .* tt .*v2 / 2 ./ R + 2 * t .* dttVi .*v2 / 2 ./ R - 2 * t .* tt .*v2 / 2 ./ R.^2 .* dRVi + 2 * t .* tt / 2 ./ R + (2*tt).*dttVi );
        end
        dtau( :, iVIndex ) = dtauVi;
    end
    
    funGradArray = dtau / b6;
end

