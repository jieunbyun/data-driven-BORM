%{
17 May 21 
Ji-Eun Byun

Sensitivity with regard to optimization results of LP + uncertain material properties
%}

clear; close all;
load ex3_optim

import ex3.*
import function.*
rng(1)

nSample = size( maxMemberForceInKNArray, 2 );
uSampleArray = randn( nMember,nSample );
yieldStrengthSampleArray_sensitvity = steelYieldStrength_knM2 + steelYieldStrength_knM2*steelYieldStrengthCov*uSampleArray;


%% Sensitivity
% c.o.v. of yield strength
optimXArray = optimXArray_m2_uncertainMaterial(:,end);

[limitFunValArray, demandLoadArray] = evalLimitFun( optimXArray, xId2MemberId, maxMemberForceInKNArray, memberForceByUnitWeight, memberLength_m, steelUnitWeight_knM3, yieldStrengthSampleArray_sensitvity );

delta = 1e-6;
cov_delta = steelYieldStrengthCov + delta;
yieldStrengthSampleArray_sensitvity_delta = steelYieldStrength_knM2 + steelYieldStrength_knM2*cov_delta*uSampleArray;
limitFunValArray_delta = evalLimitFun( optimXArray, xId2MemberId, maxMemberForceInKNArray, memberForceByUnitWeight, memberLength_m, steelUnitWeight_knM3, yieldStrengthSampleArray_sensitvity_delta );

sensitivityArray_sigma_numeric = zeros( nMember,1); sensitivityArray_sigma = zeros( nMember,1);
for iMemberIndex = 1:nMember
    iBpf = evalBpf( limitFunValArray(iMemberIndex,:) );
    
    if iBpf
        iDBpf_dSigma = limitFunGradEx3( maxMemberForceInKNArray, optimXArray, iMemberIndex, 'sigma', xId2MemberId, memberForceByUnitWeight, memberLength_m, steelUnitWeight_knM3, steelYieldStrength_knM2, yieldStrengthSampleArray_sensitvity, steelYieldStrengthCov, steelYieldStrength_knM2 );
        sensitivityArray_sigma( iMemberIndex ) = iDBpf_dSigma;

        iBpf_delta = evalBpf( limitFunValArray_delta(iMemberIndex,:) );
        iDBpf_dSigma_numeric = (iBpf_delta-iBpf)/delta;
        sensitivityArray_sigma_numeric( iMemberIndex ) = iDBpf_dSigma_numeric;
    end
end

disp( '%%% Sensitivity w.r.t. c.o.v.: Derivation, numerical, and relative difference %%%' )
disp( num2str( sensitivityArray_sigma.') )
disp( num2str( sensitivityArray_sigma_numeric.') )
disp( num2str( abs(sensitivityArray_sigma.'-sensitivityArray_sigma_numeric.')./sensitivityArray_sigma.' ) )


% Decision variable
delta = 1e-8;
sensitivityArray_x_numeric = zeros( nX, nMember ); sensitivityArray_x = zeros( nX, nMember );
for iXIndex = 1:nX
    iParam = strcat( 'x',num2str(iXIndex) );
    
    iX_delta = optimXArray; iX_delta( iXIndex ) = iX_delta( iXIndex ) + delta;
    iLimitFunValArray_delta = evalLimitFun( iX_delta, xId2MemberId, maxMemberForceInKNArray, memberForceByUnitWeight, memberLength_m, steelUnitWeight_knM3, yieldStrengthSampleArray_sensitvity );
     
    for jMemberIndex = 1:nMember
        jBpf = evalBpf( limitFunValArray(jMemberIndex,:) );

        if jBpf
            ijDBpf_dX = limitFunGradEx3( maxMemberForceInKNArray, optimXArray, jMemberIndex, iParam, xId2MemberId, memberForceByUnitWeight, memberLength_m, steelUnitWeight_knM3, steelYieldStrength_knM2, yieldStrengthSampleArray_sensitvity, steelYieldStrengthCov, steelYieldStrength_knM2, demandLoadArray );
            sensitivityArray_x( iXIndex,jMemberIndex ) = ijDBpf_dX;

            ijBpf_delta = evalBpf( iLimitFunValArray_delta(jMemberIndex,:) );
            iDBpf_dX_numeric = (ijBpf_delta-jBpf)/delta;
            sensitivityArray_x_numeric( iXIndex,jMemberIndex ) = iDBpf_dX_numeric;
        end
    end
end
sensitivityArray_x_relativeDifference = abs(sensitivityArray_x-sensitivityArray_x_numeric)./abs(sensitivityArray_x);
sensitivityArray_x_relativeDifference( abs(sensitivityArray_x) < 1e-7 ) = nan;

save ex3_sensitivity

%% Figure
fontsize_label = 16; fontName = 'times new roman'; fontsize_tick = 14;

figure('Renderer', 'painters', 'Position', [10 10 900 350])
sensitivityFigure = image( sensitivityArray_x );
sensitivityFigure.CDataMapping = 'scaled';
colorbar

xlabel( 'Member ID','Fontsize',fontsize_label,'FontName', fontName )
ylabel( 'Decision variables','Fontsize',fontsize_label+2,'FontName', fontName )
xticks(1:nMember )
yticklabels({'{\it X_1}' '{\it X_2}' '{\it X_3}' '{\it X_4}' '{\it X_5}' '{\it X_6}'})
set(gca, 'FontSize', fontsize_tick,'FontName',fontName)

saveas( gcf, 'figure/ex3_sensitivity_x.emf' )
saveas( gcf, 'figure/ex3_sensitivity_x.pdf' )


figure('Renderer', 'painters', 'Position', [10 10 900 125])
sensitivityFigure = image( sensitivityArray_sigma.' );
sensitivityFigure.CDataMapping = 'scaled';
colorbar

xlabel( 'Member ID','Fontsize',fontsize_label,'FontName', fontName )
ylabel( '{\it \delta_y}','Fontsize',fontsize_label+1,'FontName', fontName )
xticks(1:nMember )
yticklabels('')
set(gca, 'FontSize', fontsize_tick,'FontName',fontName)

yLabel = get(gca,'ylabel');
set(yLabel, 'rotation',0, 'verticalalignment', 'middle', 'HorizontalAlignment','center')
yLabel.Position(1) = yLabel.Position(1)-0.25;

saveas( gcf, 'figure/ex3_sensitivity_cov.emf' )
saveas( gcf, 'figure/ex3_sensitivity_cov.pdf' )