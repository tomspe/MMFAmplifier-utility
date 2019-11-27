% Correction: do not take conjugate of field when decomposing on mode base at each step
function [ SignalAmplitudesAlongZ ] = SignalPropagator_Corrected(ModeFieldsStack, BetasList, InjectedComplexAmplitudeList, Level2PopulationThroughoutVolume, Xaxis, Zaxis, DopingType, DopingFactor) 


switch DopingType
    case 'Erbium'
        [ ~, ~, ~, SignalGainCoefficientForTotalInversion, SignalAbsorptionCoefficientFor0Inversion ] = ErbiumAbsorptionAndGainModel(DopingFactor) ;
    case 'Ytterbium'
        [ ~, ~, ~, SignalGainCoefficientForTotalInversion, SignalAbsorptionCoefficientFor0Inversion ] = YtterbiumAbsorptionAndGainModel(DopingFactor) ;
end

NumOfModes = length(BetasList) ;
Sx = length(Xaxis) ; Sy = Sx ;
BetaStack = DuplicateValuesInListToStack(BetasList, Sy, Sx) ;  

dz = Zaxis(2) - Zaxis(1) ;
Nsteps = length(Zaxis) - 1 ;

% Allocating the output :
SignalAmplitudesAlongZ = zeros(NumOfModes, Nsteps + 1) ;
SignalAmplitudesAlongZ(:,1) = InjectedComplexAmplitudeList.' ;

%figure('color', 'w') ;
for z = 2:length(Zaxis)
    n2 = Level2PopulationThroughoutVolume(:,:,z-1) ;
    n1 = 1 - n2 ;
    gProfile = SignalGainCoefficientForTotalInversion*n2 - SignalAbsorptionCoefficientFor0Inversion*n1 ;
    FieldAmplificationProfile = (1 + gProfile*1e-6*dz/2) ;   % the factor 2 is supposed to change from power attenuation to field attenuation 
    
    % advancing the modes by dz and applying the absorption/gain profile
    AmpStack = DuplicateValuesInListToStack(SignalAmplitudesAlongZ(:,z-1), Sy, Sx) ;
    SignalField_out = FieldAmplificationProfile.*SuperposeModeFields(ModeFieldsStack, BetaStack, AmpStack, dz) ;
    % decompose this on all possible pump modes:
    Overlaps = ModeFieldsStack.*repmat(SignalField_out, 1, 1, NumOfModes) ;
    NewAmplitudes = squeeze(sum(sum(Overlaps, 1), 2)) ;
    % update output
    SignalAmplitudesAlongZ(:,z) = NewAmplitudes ; 
%     if mod(z-2,20) == 0
%         F = SuperposeModeFields(ModeFieldsStack, BetaStack, DuplicateValuesInListToStack(NewAmplitudes, Sy, Sx), 0) ;
%         imagesc(F.*conj(F)) ; colorbar ; pause(0.01) ;
%     end
    
     %disp([ num2str(100*(z - 1)/(length(Zaxis) - 1))  '% done' ])
end

end

