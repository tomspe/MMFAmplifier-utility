% Correction: do not take conjugate of field when decomposing on mode base at each step
function [ PumpAmplitudesAlongZ, Level2PopulationThroughoutVolume ] = PumpPropagator_Corrected(...
            ModeFieldsStack, BetaList, InjectedComplexAmplitudeList, Xaxis, Zaxis, dA, DopingType, DopingFactor, PlotFlag) % MultiplicationFactorToNormalizeIp, DopingType, DopingFactor, PlotFlag)
        
NumOfModes = length(BetaList) ;
Sx = length(Xaxis) ; Sy = Sx ;
BetaStack = DuplicateValuesInListToStack(BetaList, Sy, Sx) ;    
Nsteps = length(Zaxis) - 1 ;
dz = Zaxis(2) - Zaxis(1) ;
AcoreOverAtotal = pi*1^2/((2*Xaxis(end))^2) ;

[ XX, YY ] = meshgrid(Xaxis) ;
InCoreIndices = sqrt(XX.^2 + YY.^2) <= 1 ;

switch DopingType
    case 'Erbium'
        [ PumpAbsorptionCoefficient, PumpSaturationIntensity, MaximalInversionAchievable, ~, ~ ] = ErbiumAbsorptionAndGainModel(DopingFactor) ;
    case 'Ytterbium'
        [ PumpAbsorptionCoefficient, PumpSaturationIntensity, MaximalInversionAchievable, ~, ~ ] = YtterbiumAbsorptionAndGainModel(DopingFactor) ;
end
% The typical length of absorption
alpha = PumpAbsorptionCoefficient ;   %in 1/m
% The power, in theunit area in our modeling of the space (dA), which produces a saturation intensity
PumpSaturationPower = PumpSaturationIntensity*dA ;   %in W

% Allocating the output :
PumpAmplitudesAlongZ = zeros(NumOfModes, Nsteps + 1) ;
PumpAmplitudesAlongZ(:,end) = InjectedComplexAmplitudeList ;  % pump moves "backwards" - start filling from the end
Level2PopulationThroughoutVolume = zeros(Sy, Sx, Nsteps) ;
if PlotFlag
    h_fig = figure('color', 'w') ;
    h_axes = axes ;
    PlottingStep = 200 ;
    c = 1 ;
end
for z = 2:length(Zaxis)
    % Finding pump speckle pattern incoming upon this step:
    AmpStack = DuplicateValuesInListToStack(PumpAmplitudesAlongZ(:,Nsteps + 1 + 2 - z), Sy, Sx) ;  
    
        PumpField_in = SuperposeModeFields(ModeFieldsStack, BetaStack, AmpStack, 0) ;   %sum(AmpStack.*ModeFieldsStack, 3) ;
        % from field to real Intensity in W/m^2
        Fsquared = PumpField_in.*conj(PumpField_in) ;
                %PumpIntensity = MultiplicationFactorToNormalizeIp*Fsquared ; 
                % the saturation curve which translates Ip to a population of level 2
                %n2 = MaximalInversionAchievable*PumpIntensity./(PumpIntensity + PumpSaturationIntensity ) ;
        n2 = MaximalInversionAchievable*Fsquared./(Fsquared + PumpSaturationPower ) ;
        % already we can update that will be used for the signal propagator -
        % the excitation of level 2 (note that we fill this array "in reverse" - because the pump is counter-propagating) 
        Level2PopulationThroughoutVolume(:,:,Nsteps - z + 2) = n2.*InCoreIndices ;    


        % Level 2 population determines a Population Difference term
        % (((1 - r)/r)*n2 - n1) going from -1 (for no pumping) to almost 0 (saturation of pumping)
                %Pdiff = -PumpSaturationIntensity./(PumpIntensity + PumpSaturationIntensity ) ;
            %Pdiff = 2*n2 - 1 ;
        aProfile = alpha*(n2/MaximalInversionAchievable - 1)  ;
        % the population difference determines the local loss for the intensity (small absorption exp(-eps) is approximated as 1+eps):
            %FieldAttenuationProfile = (1 + alpha*1e-6*dz*Pdiff/2) ;   % the factor 2 is supposed to change from power attenuation to field attenuation
        FieldAttenuationProfile = (1 + aProfile*1e-6*dz/2) ;   % the factor 2 is supposed to change from power attenuation to field attenuation 
        % for debug: 
        %disp([ 'Mean Gain/Absorption factor: ' num2str(mean2(FieldAttenuationProfile))  '  Mean Intensity: ' num2str(mean2(PumpIntensity)*1e-9)    ])

    % Now, the field that will propagate through this step is evaluated by
    % advancing the phases dz, applying the absorption, and decomposing the
    % result upon the modes (re-coupling into the fiber)
    PumpField_out = FieldAttenuationProfile.*SuperposeModeFields(ModeFieldsStack, BetaStack, AmpStack, dz) ;   
    % decompose this on all possible pump modes:
    Overlaps = ModeFieldsStack.*repmat(PumpField_out, 1, 1, NumOfModes) ;
    NewAmplitudes = squeeze(sum(sum(Overlaps, 1), 2)) ;
    % update output that shows evolution of the pump - from the perspective of the pump itself (z not "in reverse") 
    PumpAmplitudesAlongZ(:,Nsteps + 2 - z) = NewAmplitudes ; 
     
     if PlotFlag
         if mod(z-2,PlottingStep) == 0
            imagesc(n2.*InCoreIndices,  'Parent', h_axes, [0 1]) ;  title('Population') ; pause(0.1)  ;
            MeanExcitation = (1/AcoreOverAtotal)*mean2(n2/MaximalInversionAchievable) ;
            disp([ 'Pump propagation - ' num2str(100*(z - 1)/(length(Zaxis) - 1))  '% done' ])
            c = c + 1 ;
         end
     end


end

if PlotFlag
    close(h_fig)
end

end