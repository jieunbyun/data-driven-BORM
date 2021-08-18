%{
4 May 2021
Ji-Eun Byun

Illustration of quantile, superquantile, failure probability, and buffered failure probability
%}

close all; clear;

%%
nXArray = 1e3; minX = -4; maxX = 2;
xArray = linspace( minX, maxX, nXArray );

normMean = -1; normStd = 1;
normPdf = normpdf( xArray, normMean, normStd );
normCdf = normcdf( xArray, normMean, normStd );


lineWidth = 2; fontsize_tick = 13; fontName = 'times new roman';
figure;
plot( xArray, normPdf, 'linewidth', lineWidth )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
axis( [minX maxX 0 0.45] )
saveas( gcf,  'figure/normPdf.emf' )
saveas( gcf, 'figure/normPdf.pdf' )

figure;
plot( xArray, normCdf, 'linewidth', lineWidth )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
saveas( gcf, 'figure/normCdf.emf' )
saveas( gcf, 'figure/normCdf.pdf' )
ax.FontName = fontName;
axis( [minX maxX 0 1] )


truncatedMeanArray = normMean + normStd * normpdf( ( xArray - normMean )/normStd, 0, 1 ) ./ ( 1 - normcdf( ( xArray - normMean )/normStd, 0, 1 ) );
figure;
plot( xArray, truncatedMeanArray, 'linewidth', lineWidth )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
saveas( gcf, 'figure/truncatedMean.emf' )
saveas( gcf, 'figure/truncatedMean.pdf' )
ax.FontName = fontName;
% axis( [minX maxX 0 1] )

%% Quantile, superquantile, failure probability, buffered failure probability
quantile0 = 0;
pf = 1- normcdf( quantile0, normMean, normStd );
superquantile0 = normMean + normStd * normpdf( (quantile0-normMean)/normStd ) / ( 1-pf );
bpf = 1 - normcdf( superquantile0, normMean, normStd );
disp( ['pf: ' num2str( pf )] )
disp( ['bpf: ' num2str( bpf )] )
disp( ['quantile: ' num2str( quantile0 )] )
disp( ['superquantile: ' num2str( superquantile0 )] )