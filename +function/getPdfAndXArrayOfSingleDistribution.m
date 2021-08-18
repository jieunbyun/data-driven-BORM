function [pdfArrayOfSingleDistribution, xArrayOfSingleDistribution] = getPdfAndXArrayOfSingleDistribution( pf, pdfFun, rndFun, xArray, nSample )

pdfArrayOfSingleDistribution = pdfFun( xArray );

sampleArray = rndFun( nSample );
sampleArraySorted = sort( sampleArray, 'descend' );

thresholdSampleIndex = ceil( nSample * pf );
thresholdSampleValue = sampleArraySorted( thresholdSampleIndex );
xArrayOfSingleDistribution = xArray - thresholdSampleValue;

