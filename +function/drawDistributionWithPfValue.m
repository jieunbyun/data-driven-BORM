function pdfArray = drawDistributionWithPfValue( pf, xArray, paramArray, distributionIndex, linspec, lineWidthArray, markerSizeArray, colorSpec, fontsize_tick, fontName, fontsize_label, axisRange, legendParam )

import function.*

mu = 1; % arbitrary positive number
nSample = 1e5;

xArray = xArray(:).';
nX = length( xArray );
nDist = length( paramArray );
pdfArray = zeros( nDist, nX );
xArrayOfDistributions = zeros( nDist, nX );


normPdfFun = @(xArray) normpdf( xArray, mu, paramArray( distributionIndex.normal ) );
normRndFun = @( nSample ) normrnd(  mu, paramArray( distributionIndex.normal ), nSample, 1 );
[pdfArray( distributionIndex.normal, : ), xArrayOfDistributions( distributionIndex.normal, : )] = getPdfAndXArrayOfSingleDistribution( pf, normPdfFun, normRndFun, xArray, nSample );

expPdfFun = @(xArray) exppdf( xArray, 1/mu );
expRndFun = @( nSample ) exprnd( 1/mu, nSample, 1 );
[pdfArray( distributionIndex.exp, : ), xArrayOfDistributions( distributionIndex.exp, : )] = getPdfAndXArrayOfSingleDistribution( pf, expPdfFun, expRndFun, xArray, nSample );

paramLogn = paramArray( distributionIndex.logn );
lognPdfFun = @(xArray) lognpdf( xArray, -paramLogn^2/2, paramLogn );
lognRndFun = @(nSample) lognrnd( -paramLogn^2/2, paramLogn, nSample, 1 );
[pdfArray( distributionIndex.logn, : ), xArrayOfDistributions( distributionIndex.logn, : )] = getPdfAndXArrayOfSingleDistribution( pf, lognPdfFun, lognRndFun, xArray, nSample );

paramWeibull = paramArray( distributionIndex.weibull );
weibullPdfFun = @(xArray) wblpdf( xArray, 1/gamma( 1+1/paramWeibull ), paramWeibull );
weibullRndFun = @(nSample) wblrnd( 1/gamma( 1+1/paramWeibull ), paramWeibull, nSample, 1 );
[pdfArray( distributionIndex.weibull, : ), xArrayOfDistributions( distributionIndex.weibull, : )] = getPdfAndXArrayOfSingleDistribution( pf, weibullPdfFun, weibullRndFun, xArray, nSample );

paramGev = paramArray( distributionIndex.gev ); sigma_gev = 1; mu_gev = 1-sigma_gev*double( eulergamma );
gevPdfFun = @(xArray) gevpdf( xArray, paramGev, sigma_gev, mu_gev );
gevRndFun = @(nSample) gevrnd( paramGev, sigma_gev, mu_gev, nSample, 1);
[pdfArray( distributionIndex.gev, : ), xArrayOfDistributions( distributionIndex.gev, : )] = getPdfAndXArrayOfSingleDistribution( pf, gevPdfFun, gevRndFun, xArray, nSample );

figure;
grid on
hold on
for iDistIndex = 1:nDist
    plot( xArrayOfDistributions( iDistIndex,: ), pdfArray( iDistIndex,: ), linspec{iDistIndex}, 'linewidth', lineWidthArray( iDistIndex ), 'markersize', markerSizeArray( iDistIndex ), 'color', colorSpec( iDistIndex,: ) )
end
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
xlabel( '\it{y}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it f_Y} ({\it y})', 'fontname', fontName, 'fontsize', fontsize_label )
axis( axisRange.wholeDistribution )
legend( legendParam.labels, 'fontname', legendParam.font, 'fontsize', legendParam.fontSize, 'location', legendParam.location )


figure;
grid on
hold on
for iDistIndex = 1:nDist
    plot( xArrayOfDistributions( iDistIndex,: ), pdfArray( iDistIndex,: ), linspec{iDistIndex}, 'linewidth', lineWidthArray( iDistIndex ), 'markersize', markerSizeArray( iDistIndex ), 'color', colorSpec( iDistIndex,: ) )
end
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
xlabel( '\it{y}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it f_Y} ({\it y})', 'fontname', fontName, 'fontsize', fontsize_label )
axis( axisRange.rightTail )
legend( legendParam.labels, 'fontname', legendParam.font, 'fontsize', legendParam.fontSize, 'location', legendParam.location )