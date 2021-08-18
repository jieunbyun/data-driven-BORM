function [pfArray, tauArray] = evalPfAndTauFromBpfAndSample( bpfArray, sampleArray )

    pfArray = zeros( size( bpfArray ) );
    for iPfIndex = 1:length( pfArray )
        iBpf = bpfArray( iPfIndex );
        iPf = evalPf1FromBpf1AndSample( iBpf, sampleArray );
        pfArray( iPfIndex ) = iPf;
    end
    
    tauArray = bpfArray ./ pfArray;
    
end

function pf1 = evalPf1FromBpf1AndSample( bpf1, sampleArray )

   z0 = quantile( sampleArray, 1-bpf1 );
   zArray = sampleArray - z0; zArray = zArray( zArray > 0 );
    
    alphaQuantile = mean( zArray ) + z0;
    pf1 = mean( sampleArray > alphaQuantile );

end