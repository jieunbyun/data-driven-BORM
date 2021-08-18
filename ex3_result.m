%{
13 May 21 
Ji-Eun Byun

Results of optimization
%}

clear; close all;
load ex3_optim

import ex3.*
import function.*
rng(1)

nSample = size( maxMemberForceInKNArray, 2 );

%% LP
optimXArray_LP = optimXArray_m2(:,end);
limitFunValArray_LP = evalLimitFun( optimXArray_LP, xId2MemberId, maxMemberForceInKNArray, memberForceByUnitWeight, memberLength_m, steelUnitWeight_knM3, steelYieldStrength_knM2 );

bPfArray_LP = zeros( 1,nMember ); pfArray_LP = zeros( 1,nMember );
for iMemberIndex = 1:nMember
    iLimitFunValArray = limitFunValArray_LP( iMemberIndex,: );
    
    pfArray_LP( iMemberIndex ) = mean( iLimitFunValArray > 0 );
    bPfArray_LP( iMemberIndex ) = evalBpf( iLimitFunValArray );
end
tauArray_LP = bPfArray_LP./pfArray_LP;
tauArray_LP( isnan(tauArray_LP ) ) = 0;

disp( '%%%%% Results of LP  %%%%%' )
disp( ['Bpf: ' num2str(bPfArray_LP)] )
disp( ['Pf: ' num2str(pfArray_LP)] )
disp( ['tau: ' num2str(bPfArray_LP./pfArray_LP)] )


%% LP + uncertain material properties
uSampleArray = normrnd( 0, 1, [nMember, nSample] );
steelYieldStrengthSampleArray = steelYieldStrength_knM2 + steelYieldStrength_knM2 * steelYieldStrengthCov * uSampleArray;

optimXArray_LPAndUncertainMaterial = optimXArray_m2_uncertainMaterial(:,end);
limitFunValArray_LPAndUncertainMaterial = evalLimitFun( optimXArray_LPAndUncertainMaterial, xId2MemberId, maxMemberForceInKNArray, memberForceByUnitWeight, memberLength_m, steelUnitWeight_knM3, steelYieldStrengthSampleArray );

bPfArray_LPAndUncertainMaterial = zeros( 1,nMember ); pfArray_LPAndUncertainMaterial = zeros( 1,nMember );
for iMemberIndex = 1:nMember
    iLimitFunValArray = limitFunValArray_LPAndUncertainMaterial( iMemberIndex,: );
    
    pfArray_LPAndUncertainMaterial( iMemberIndex ) = mean( iLimitFunValArray > 0 );
    bPfArray_LPAndUncertainMaterial( iMemberIndex ) = evalBpf( iLimitFunValArray );
end
tauArray_LPAndUncertainMaterial = bPfArray_LPAndUncertainMaterial./pfArray_LPAndUncertainMaterial;
tauArray_LPAndUncertainMaterial( isnan(tauArray_LPAndUncertainMaterial ) ) = 0;

disp( '%%%%% Results of LP + uncertain material property %%%%%' )
disp( ['Bpf: ' num2str(bPfArray_LPAndUncertainMaterial)] )
disp( ['Pf: ' num2str(pfArray_LPAndUncertainMaterial)] )
disp( ['tau: ' num2str(bPfArray_LPAndUncertainMaterial./pfArray_LPAndUncertainMaterial)] )

%% MIP + uncertainty in material property
optimXArray_MIP = optimXArray_m2_uncertainMaterial_milp(:,end);
limitFunValArray_MIP = evalLimitFun( optimXArray_MIP, xId2MemberId, maxMemberForceInKNArray, memberForceByUnitWeight, memberLength_m, steelUnitWeight_knM3, steelYieldStrengthSampleArray );

bPfArray_MIP = zeros( 1,nMember ); pfArray_MIP = zeros( 1,nMember );
for iMemberIndex = 1:nMember
    iLimitFunValArray = limitFunValArray_MIP( iMemberIndex,: );
    
    pfArray_MIP( iMemberIndex ) = mean( iLimitFunValArray > 0 );
    bPfArray_MIP( iMemberIndex ) = evalBpf( iLimitFunValArray );
end
tauArray_MIP = bPfArray_MIP./pfArray_MIP;
tauArray_MIP( isnan(tauArray_MIP ) ) = 0;

disp( '%%%%% Results of MIP %%%%%' )
disp( ['Bpf: ' num2str(bPfArray_MIP)] )
disp( ['Pf: ' num2str(pfArray_MIP)] )
disp( ['tau: ' num2str(bPfArray_MIP./pfArray_MIP)] )

save ex3_result

%% Figure
markerSize = 7; fontsize_tick = 14; fontName = 'times new roman'; fontsize_label = 16; lineWidth = 1.5; fontsize_legend = 14;
figure('Renderer', 'painters', 'Position', [10 10 900 350])
plot( 1:nMember, bPfArray_LP, 's', 'linewidth', lineWidth, 'markerSize', markerSize+1 )
hold on
plot( 1:nMember, bPfArray_LPAndUncertainMaterial, 'o', 'linewidth', lineWidth, 'markerSize', markerSize-2 )
plot( 1:nMember, bPfArray_MIP, '*', 'linewidth', lineWidth-0.5, 'markerSize', markerSize )
% set( gca, 'YScale', 'log' )
% axis( [12.58 12.72 3e-3 1] )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
ax.XTick = 1:nMember;
xlabel( 'Member ID', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( 'Buffered failure probability', 'fontname', fontName, 'fontsize', fontsize_label )
legend( {'Continuous' 'Continuous & Uncertain {\it \sigma_y}', 'Discrete & Uncertain {\it \sigma_y}'}, 'fontname', fontName, 'fontsize', fontsize_legend )

saveas( gcf, 'figure/ex3_bpf.emf' )
saveas( gcf, 'figure/ex3_bpf.pdf' )


figure('Renderer', 'painters', 'Position', [10 10 900 350])
plot( 1:nMember, pfArray_LP, 's', 'linewidth', lineWidth, 'markerSize', markerSize+1 )
hold on
plot( 1:nMember, pfArray_LPAndUncertainMaterial, 'o', 'linewidth', lineWidth, 'markerSize', markerSize-2 )
plot( 1:nMember, pfArray_MIP, '*', 'linewidth', lineWidth-0.5, 'markerSize', markerSize )
% set( gca, 'YScale', 'log' )
% axis( [12.58 12.72 3e-3 1] )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
ax.XTick = 1:nMember;
xlabel( 'Member ID', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( 'Conventional failure probability', 'fontname', fontName, 'fontsize', fontsize_label )

saveas( gcf, 'figure/ex3_pf.emf' )
saveas( gcf, 'figure/ex3_pf.pdf' )


figure('Renderer', 'painters', 'Position', [10 10 900 350])
plot( 1:nMember, tauArray_LP, 's', 'linewidth', lineWidth, 'markerSize', markerSize+1 )
hold on
plot( 1:nMember, tauArray_LPAndUncertainMaterial, 'o', 'linewidth', lineWidth, 'markerSize', markerSize-2 )
plot( 1:nMember, tauArray_MIP, '*', 'linewidth', lineWidth-0.5, 'markerSize', markerSize )
% set( gca, 'YScale', 'log' )
% axis( [12.58 12.72 3e-3 1] )
grid on
ax = gca;
ax.FontSize = fontsize_tick;
ax.FontName = fontName;
ax.XTick = 1:nMember;
xlabel( 'Member ID', 'fontname', fontName, 'fontsize', fontsize_label )
ylabel( 'Buffered tail index', 'fontname', fontName, 'fontsize', fontsize_label )

saveas( gcf, 'figure/ex3_tau.emf' )
saveas( gcf, 'figure/ex3_tau.pdf' )