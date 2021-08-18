%{
6 May 21
Ji-Eun Byun

Compute buffered tail index (tau) for common distributions
%}

close all; clear
import function.*
rng(1)

lineWidth = 2; fontsize_tick = 12; fontName = 'times new roman'; fontsize_label = 16; fontsize_legend = 13;

bpfArray = 10.^( -5:0.1:0 );
nBpf = length( bpfArray );
alphaArray = 1 - bpfArray;
xAxisRange = [1e-6 1];

nSample = 1e7;
bpfArray_sample = 10.^( -4:0.1:0 );
nBpf_sample = length( bpfArray_sample );

%% Normal distribution
% Analytical derivation
pfArray_normal = ( 1 - normcdf( normpdf( sqrt(2)*erfinv( 2*alphaArray - 1 ) ) ./ bpfArray ) );
tau_normal = bpfArray ./ pfArray_normal;

% By sampling
sampleArray = normrnd( 0, 1, [nSample, 1] );
[pfArray_normal_sample, tau_normal_sample] = evalPfAndTauFromBpfAndSample( bpfArray_sample, sampleArray );

figure;
plot( pfArray_normal, tau_normal, 'linewidth', lineWidth )
set(gca, 'XScale', 'log');
hold on
plot( pfArray_normal_sample, tau_normal_sample, ':', 'linewidth', lineWidth+1 )
plot( [1e-6 0.01 0.3 0.5], [2.68 2.61 2.4 2], ':', 'linewidth', lineWidth )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
axis( [xAxisRange 2 2.7] )
xlabel( '\it{p_f}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it \tau}_{normal}', 'fontname', fontName, 'fontsize', fontsize_label )
legend( {'Analytical derivation' '10^7 samples' 'Proposed function of \tau'}, 'fontname', fontName, 'fontsize', fontsize_legend, 'location', 'best' )
% saveas( gcf,  'figure/tau_normal.emf' )
% saveas( gcf, 'figure/tau_normal.pdf' )


%% Exponential distribution
% Analytical derivation
tau_exp = exp(1) * ones( size( bpfArray ) );

% By sampling
sampleArray = exprnd( 1, nSample, 1 );
[pfArray_exp_sample, tau_exp_sample] = evalPfAndTauFromBpfAndSample( bpfArray_sample, sampleArray );

figure;
plot( bpfArray./tau_exp, tau_exp, 'linewidth', lineWidth )
set(gca, 'XScale', 'log');
hold on
plot( pfArray_exp_sample, tau_exp_sample, ':', 'linewidth', lineWidth+1 )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
axis( [xAxisRange 2.5 3] )
xlabel( '\it{p_f}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it \tau}_{exp}', 'fontname', fontName, 'fontsize', fontsize_label )
legend( {'Derivation' 'Sampling'}, 'fontname', fontName, 'fontsize', fontsize_legend, 'location', 'best' )

%% Lognormal distribution
logStd = 0.125;
pfArray_lognormal = 1 - normcdf( -log(2)/logStd + logStd/2 + 1/logStd*log( (1+erf( logStd/sqrt(2) - erfinv( 2*alphaArray - 1 ) ) ) ./ (1-alphaArray) ) );
tau_lognormal = bpfArray ./ pfArray_lognormal;

% By sampling
sampleArray = lognrnd( 0, logStd, [nSample, 1] );
[pfArray_lognormal_sample, tau_lognormal_sample] = evalPfAndTauFromBpfAndSample( bpfArray_sample, sampleArray );

figure;
plot( pfArray_lognormal, tau_lognormal, 'linewidth', lineWidth )
set(gca, 'XScale', 'log');
hold on
plot( pfArray_lognormal_sample, tau_lognormal_sample, ':', 'linewidth', lineWidth+1 )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
xlabel( '\it{p_f}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it \tau}_{lognormal}', 'fontname', fontName, 'fontsize', fontsize_label )
legend( {'Derivation' 'Sampling'}, 'fontname', fontName, 'fontsize', fontsize_legend, 'location', 'best' )

% % Different standard deviations
logStdArray = [0.125 1];
nLogStd = length( logStdArray );
pfArray_lognormal_param = zeros( nLogStd, nBpf ); tau_lognormal_param = zeros( nLogStd, nBpf );
figure;
hold on
for iStdIndex = 1:nLogStd
    iLogStd = logStdArray( iStdIndex );
    iPfArray_lognormal = 1 - normcdf( -log(2)/iLogStd + iLogStd/2 + 1/iLogStd*log( (1+erf( iLogStd/sqrt(2) - erfinv( 2*alphaArray - 1 ) ) ) ./ (1-alphaArray) ) );
    iTau_lognormal = bpfArray ./ iPfArray_lognormal;
    
    plot( iPfArray_lognormal, iTau_lognormal, 'linewidth', lineWidth )
    
    pfArray_lognormal_param( iStdIndex, : ) = iPfArray_lognormal;
    tau_lognormal_param( iStdIndex, : ) = iTau_lognormal;
end
set(gca, 'XScale', 'log');
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
xlabel( '\it{p_f}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it \tau}_{lognormal}', 'fontname', fontName, 'fontsize', fontsize_label )

%% Weibull
shapeParameter = 0.5;
pfArray_weibull = exp( - (1./bpfArray .* gammainc( -log( bpfArray ), 1+1/shapeParameter, 'upper' ) .* gamma( 1+1/shapeParameter ) ).^shapeParameter );
tau_weibull = bpfArray ./ pfArray_weibull;

sampleArray = wblrnd( 1, shapeParameter, [nSample, 1] );
[pfArray_weibull_sample, tau_weibull_sample] = evalPfAndTauFromBpfAndSample( bpfArray_sample, sampleArray );

figure;
plot( pfArray_weibull, tau_weibull, 'linewidth', lineWidth )
set(gca, 'XScale', 'log');
hold on
plot( pfArray_weibull_sample, tau_weibull_sample, ':', 'linewidth', lineWidth+1 )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
xlabel( '\it{p_f}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it \tau}_{Weibull}', 'fontname', fontName, 'fontsize', fontsize_label )
legend( {'Derivation' 'Sampling'}, 'fontname', fontName, 'fontsize', fontsize_legend, 'location', 'best' )

% % Different parameters
shapeParameterArray = [0.5 1.5];
nShapeParameter = length( shapeParameterArray );
pfArray_weibull_param = zeros( nShapeParameter, nBpf ); tau_weibull_param = zeros( nShapeParameter, nBpf );
figure;
hold on
for iParamIndex = 1:nShapeParameter
    iShapeParam = shapeParameterArray( iParamIndex );
    iPfArray = exp( - (1./bpfArray .* gammainc( -log( bpfArray ), 1+1/iShapeParam, 'upper' ) .* gamma( 1+1/iShapeParam ) ).^iShapeParam );
    iTau = bpfArray ./ iPfArray;
    
    plot( iPfArray, iTau, 'linewidth', lineWidth )
    
    pfArray_weibull_param( iParamIndex, : ) = iPfArray;
    tau_weibull_param( iParamIndex, : ) = iTau;
end
set(gca, 'XScale', 'log');
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
xlabel( '\it{p_f}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it \tau}_{Weibull}', 'fontname', fontName, 'fontsize', fontsize_label )

%% GEV
epsilon = 0;
if ~epsilon
    pfArray_gev = 1 - exp( -exp( -1./bpfArray.*( double( eulergamma ) + alphaArray .* log( -log( alphaArray ) ) - logint( alphaArray ) ) ) );
else
    pfArray_gev = 1 - exp( -( 1 + 1./bpfArray .* ( gammainc( log( 1./alphaArray ), 1-epsilon ) * gamma( 1-epsilon ) - bpfArray ) ) .^ ( -1/epsilon ) );
end
tau_gev = bpfArray ./ pfArray_gev;

sampleArray = gevrnd( epsilon, 1, 0, [nSample, 1] );
[pfArray_gev_sample, tau_gev_sample] = evalPfAndTauFromBpfAndSample( bpfArray_sample, sampleArray );

figure;
plot( pfArray_gev, tau_gev, 'linewidth', lineWidth )
set(gca, 'XScale', 'log');
hold on
plot( pfArray_gev_sample, tau_gev_sample, ':', 'linewidth', lineWidth+1 )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
xlabel( '\it{p_f}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it \tau}_{GEV}', 'fontname', fontName, 'fontsize', fontsize_label )
legend( {'Derivation' 'Sampling'}, 'fontname', fontName, 'fontsize', fontsize_legend, 'location', 'best' )

% % Different parameters
shapeParameterArray_GEV = [0 0.5];
nParameter_GEV = length( shapeParameterArray_GEV );
pfArray_gev_param = zeros( nParameter_GEV, nBpf ); tau_gev_param = zeros( nParameter_GEV, nBpf );
figure;
hold on
for iParamIndex = 1:nParameter_GEV
    iShapeParam = shapeParameterArray_GEV( iParamIndex );
    
    if ~iShapeParam
        iPfArray = 1 - exp( -exp( -1./bpfArray.*( double( eulergamma ) + alphaArray .* log( -log( alphaArray ) ) - logint( alphaArray ) ) ) );
    else
        iPfArray = 1 - exp( -( 1 + 1./bpfArray .* ( gammainc( log( 1./alphaArray ), 1-iShapeParam ) * gamma( 1-iShapeParam ) - bpfArray ) ) .^ ( -1/iShapeParam ) );
    end
    iTau = bpfArray ./ iPfArray;
    
    plot( iPfArray, iTau, 'linewidth', lineWidth )
    
    pfArray_gev_param( iParamIndex, : ) = iPfArray;
    tau_gev_param( iParamIndex, : ) = iTau;
end
set(gca, 'XScale', 'log');
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
xlabel( '\it{p_f}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it \tau}_{GEV}', 'fontname', fontName, 'fontsize', fontsize_label )

save demo2_bufferedTailIndex

%% Comparison of distributions
colorSpec = get(gca, 'colororder' );
% linspec = {'-', '--', 'o', 'x', 's'};
linspec = {'-', '+-', '--', ':', '-.'};
lineWidthArray = [1 0.7 1.5 2 1.7];
markerSizeArray = [1 3 4 5 5];

figure;
figIndex = 0;
figIndex = figIndex + 1;
plot( pfArray_normal, tau_normal, linspec{figIndex}, 'linewidth', lineWidthArray( figIndex ), 'color', colorSpec( figIndex,: ) )
set(gca, 'XScale', 'log');
grid on
hold on

figIndex = figIndex + 1;
plot( bpfArray./tau_exp, tau_exp, linspec{figIndex}, 'linewidth', lineWidthArray( figIndex ), 'markersize', markerSizeArray(figIndex), 'color', colorSpec( figIndex,: ) )

figIndex = figIndex + 1; 
plot( pfArray_lognormal_param(1,:), tau_lognormal_param(1,:),  linspec{figIndex}, 'linewidth', lineWidthArray( figIndex ), 'markersize', markerSizeArray(figIndex), 'color', colorSpec( figIndex,: ) )

figIndex = figIndex + 1; 
plot( pfArray_weibull_param(2,:), tau_weibull_param(2,:),  linspec{figIndex}, 'linewidth', lineWidthArray( figIndex ), 'markersize', markerSizeArray(figIndex), 'color', colorSpec( figIndex,: )  )

figIndex = figIndex + 1; 
plot( pfArray_gev_param(1,1:(end-1)), tau_gev_param(1,1:(end-1)), linspec{figIndex}, 'linewidth', lineWidthArray( figIndex ), 'markersize', markerSizeArray(figIndex), 'color', colorSpec( figIndex,: )  )

ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
xlabel( '\it{p_f}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it \tau}', 'fontname', fontName, 'fontsize', fontsize_label )
axis( [1e-6 1e0 2.2 2.8] )

legend( {'Normal', 'Exponential', 'Lognormal, {\it s} = 0.125', 'Weibull, {\it k} = 1.5', 'GEV, {\it \xi} = 0' }, 'fontname', fontName, 'fontsize', fontsize_legend-1, 'location', 'southwest' )
% saveas( gcf,  'figure/tau_nonheavy_tailed.emf' )
% saveas( gcf, 'figure/tau_nonheavy_tailed.pdf' )


figure;
figIndex = 0;
figIndex = figIndex + 1;
plot( pfArray_normal, tau_normal, linspec{figIndex}, 'linewidth', lineWidthArray( figIndex ), 'color', colorSpec( figIndex,: ) )
set(gca, 'XScale', 'log');
grid on
hold on

figIndex = figIndex + 1;
plot( bpfArray./tau_exp, tau_exp, linspec{figIndex}, 'linewidth', lineWidthArray( figIndex ), 'markersize', markerSizeArray(figIndex), 'color', colorSpec( figIndex,: ) )

figIndex = figIndex + 1; 
plot( pfArray_lognormal_param(2,:), tau_lognormal_param(2,:),  linspec{figIndex}, 'linewidth', lineWidthArray( figIndex ), 'markersize', markerSizeArray(figIndex), 'color', colorSpec( figIndex,: ) )

figIndex = figIndex + 1; 
plot( pfArray_weibull_param(1,:), tau_weibull_param(1,:), linspec{figIndex}, 'linewidth', lineWidthArray( figIndex ), 'markersize', markerSizeArray(figIndex), 'color', colorSpec( figIndex,: )  )

figIndex = figIndex + 1; 
plot( pfArray_gev_param(2,:), tau_gev_param(2,:), linspec{figIndex}, 'linewidth', lineWidthArray( figIndex ), 'markersize', markerSizeArray(figIndex), 'color', colorSpec( figIndex,: )  )

ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
xlabel( '\it{p_f}', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( '{\it \tau}', 'fontname', fontName, 'fontsize', fontsize_label )

legend( {'Normal', 'Exponential', 'Lognormal, {\it s} = 1', 'Weibull, {\it k} = 0.5', 'GEV, {\it \xi} = 0.5' }, 'fontname', fontName, 'fontsize', fontsize_legend-1, 'location', 'best' )
% saveas( gcf,  'figure/tau_heavy_tailed.emf' )
% saveas( gcf, 'figure/tau_heavy_tailed.pdf' )

%% Draw distributions
distributionIndex.normal = 1; distributionIndex.exp = 2; distributionIndex.logn = 3; distributionIndex.weibull = 4; distributionIndex.gev = 5;
xArrayNonheavyTail = -2:0.05:7; xArray_upperTail = 3:0.1:6;

pf = 1e-1;

axisRangeNonheavyTail.wholeDistribution = [-2 2 0 1.2]; axisRangeNonheavyTail.rightTail = [0 3 0 0.1];
legendParamNonheavyTail.labels = {'Normal', 'Exponential', 'Lognormal, {\it s} = 0.125', 'Weibull, {\it k} = 1.5', 'GEV, {\it \xi} = 0' }; legendParamNonheavyTail.font = fontName; legendParamNonheavyTail.fontSize = fontsize_legend-1; legendParamNonheavyTail.location = 'best';
drawDistributionWithPfValue( pf, xArrayNonheavyTail, [1 1 logStdArray(1) shapeParameterArray(2) shapeParameterArray_GEV(1)], distributionIndex, linspec, lineWidthArray, markerSizeArray, colorSpec, fontsize_tick, fontName, fontsize_label, axisRangeNonheavyTail, legendParamNonheavyTail );

% saveas( gcf,  'figure/distribution_nonheavy_tailed.emf' )
% saveas( gcf, 'figure/distribution_nonheavy_tailed.pdf' )

% saveas( gcf,  'figure/distribution_nonheavy_tailed_upperTail.emf' )
% saveas( gcf, 'figure/distribution_nonheavy_tailed_upperTail.pdf' )


xArrayHeavyTail = -10:0.1:25;
axisRangeHeavyTail.wholeDistribution = [-2.3 6 0 1.2]; axisRangeHeavyTail.rightTail = [0 15 0 0.01];
legendParamHeavyTail = legendParamNonheavyTail;
legendParamHeavyTail.labels = {'Normal', 'Exponential', 'Lognormal, {\it s} = 1', 'Weibull, {\it k} = 0.5', 'GEV, {\it \xi} = 0.5' };
drawDistributionWithPfValue( pf, xArrayHeavyTail, [1 1 logStdArray(2) shapeParameterArray(1) shapeParameterArray_GEV(2)], distributionIndex, linspec, lineWidthArray, markerSizeArray, colorSpec, fontsize_tick, fontName, fontsize_label, axisRangeHeavyTail, legendParamHeavyTail );

% saveas( gcf,  'figure/distribution_heavy_tailed.emf' )
% saveas( gcf, 'figure/distribution_heavy_tailed.pdf' )

% saveas( gcf,  'figure/distribution_heavy_tailed_upperTail.emf' )
% saveas( gcf, 'figure/distribution_heavy_tailed_upperTail.pdf' )