% Correction: do not take conjugate of field when decomposing on mode base at each step
% Also: all gain-induced contributions from some specific mode n
% should not depend on Beta(n), as was previousely assumed in the model
% based on taylor expansion of the e^j*Beta(n)z. Therefore 
function [ TransmissionMatrix, TMsAlongZ ] = CalculateTransmissionMatrix_Corrected(...
            ModeCouplingArray_s, BetaList_s, Level2PopulationThroughoutVolume, Xaxis, Zaxis, DopingType, DopingFactor, PlotFlag)


NumOfSignalModes = length(BetaList_s) ;
NumOfCouples = size(ModeCouplingArray_s, 3) ;   % should be equal to (NumOfSignalModes + 1)*NumOfSignalModes/2 
Sx = length(Xaxis) ; Sy = Sx ;
%BetaStack = DuplicateValuesInListToStack(BetaList_s, Sy, Sx) ;    
Nsteps = length(Zaxis) - 1 ;
dz = Zaxis(2) - Zaxis(1) ; 


switch DopingType
    case 'Erbium'
        [ ~, ~, MaximalInversionAchievable, SignalGainCoefficientForTotalInversion, SignalAbsorptionCoefficientFor0Inversion ] = ErbiumAbsorptionAndGainModel(DopingFactor) ;
    case 'Ytterbium'
        [ ~, ~, MaximalInversionAchievable, SignalGainCoefficientForTotalInversion, SignalAbsorptionCoefficientFor0Inversion ] = YtterbiumAbsorptionAndGainModel(DopingFactor) ;
end

% At every dz step, we will add the matrix whose elements represent the
% mode-to-mode induced gain, to this diagonal matrix which represents the
% simple propagation within the fiber (passive)
PassiveTM = diag(exp(-1i*dz*1e-6*BetaList_s)) ;
% matrix whose elements represent the mode-to-mode - we constantly have in
% each column the phaseof the contributing mode; this will be later
% multiplied, at each dz step, but an overlap factor that depends on the pump speckle for that step 
UniformGainTM = repmat(exp(-1i*dz*1e-6*BetaList_s), NumOfSignalModes, 1) ;       
            %UniformGainTM = repmat(-1i*dz*1e-6*BetaList_s.*exp(-1i*dz*1e-6*BetaList_s), NumOfSignalModes, 1) ;

% allocating the output
TransmissionMatrix = eye(size(PassiveTM)) ;
if PlotFlag
%     h_fig = figure('color', 'w') ;
%     h_axes = axes ;
    SampleEachXsteps = 20 ;
    %TMsAlongZ = zeros(NumOfSignalModes*NumOfSignalModes, ceil(Nsteps/SampleEachXsteps)) ;
    SampleCounter = 0 ;
else
    %TMsAlongZ = [] ;
end
% compare to a hypothetical spatially-uniform pump 
[ XX, YY ] = meshgrid(Xaxis) ;  InCoreIndices = sqrt(XX.^2 + YY.^2) <= 1 ;   
% n2 = MaximalInversionAchievable*InCoreIndices ;    n1 = (1 - n2).*InCoreIndices ;
% gProfile = SignalGainCoefficientForTotalInversion*n2 - SignalAbsorptionCoefficientFor0Inversion*n1 ;
% dComplexRefractiveIndexProfile = gProfile*1e-6*dz/2 ;  % the factor 2 is supposed to change from power attenuation to field attenuation 
% IntegralsArray = ModeCouplingArray_s.*repmat(dComplexRefractiveIndexProfile, 1, 1, NumOfCouples) ;
% GmatrixSUP = ArrangeListAsSymmetricMatrix(squeeze(sum(sum(IntegralsArray, 1), 2)), NumOfSignalModes) ;
% GmatrixSUP = diag(GmatrixSUP) ;

% propagation of real pump
TMsAlongZ = zeros(NumOfSignalModes, floor((length(Zaxis)-1)/2)) ;
for z = 2:length(Zaxis)
    % using what was previously calculated for the pump:
    n2 = Level2PopulationThroughoutVolume(:,:,z-1) ;
    n1 = (1 - n2).*InCoreIndices ;
    % gain/absorption is still in 1/m units
    gProfile = SignalGainCoefficientForTotalInversion*n2 - SignalAbsorptionCoefficientFor0Inversion*n1 ;
    % translating to a unitless complex n: 
    dComplexRefractiveIndexProfile = gProfile*1e-6*dz/2 ;  % the factor 2 is supposed to change from power attenuation to field attenuation 
     
    % the same n profile (induced by the pump distribution) appears identicaly in all integrals that
    % represent the mode coupling terms
    IntegralsArray = ModeCouplingArray_s.*repmat(dComplexRefractiveIndexProfile, 1, 1, NumOfCouples) ;
    % Arranging this (N+1)*N/2 vector into a NxN matrix, we get something
    % that represents the overlap integrals with a "complex refractive
    % index" factor, we still need to multiply this by the "idealized gain"
    % matrix mentioned above, and add to the passive propagation that does
    % not depend on the gain (and is therefore diagonal in our mode basis)
    Gmatrix = ArrangeListAsSymmetricMatrix(squeeze(sum(sum(IntegralsArray, 1), 2)), NumOfSignalModes) ;
    tm = PassiveTM + Gmatrix.*UniformGainTM ;
    % lastly, accumulating the effect of all slices by multiplication:
    TransmissionMatrix = tm*TransmissionMatrix ;
    if mod(z,2)
        TMsAlongZ(:,floor(z/2)) = diag(Gmatrix) ; %./GmatrixSUP ;
    end
      if PlotFlag && mod(z-2,SampleEachXsteps) == 0
          SampleCounter = SampleCounter + 1 ;
          %TMsAlongZ(:, SampleCounter) = TransmissionMatrix(:) ;
%         imagesc(abs(Gmatrix.*UniformGainTM), 'Parent', h_axes) ; colorbar ; pause(0.01) ;
           disp([ 'Signal propagation - ' num2str(100*(z - 1)/(length(Zaxis) - 1))  '% done' ])
      end
      
end


% if PlotFlag
% %    close(h_fig)
%     TMsAlongZ = TMsAlongZ(:, 1:SampleCounter) ;
% end


end