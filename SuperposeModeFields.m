

function SummedField = SuperposeModeFields(ModeFieldsStack, BetaStack, ComplexAmplitudeStack, z_in_um)

% For a given propagation distance from the fiber entrance, sum all the
% modes with the correct phase for each, to receive the total field in the
% fiber at that z point
% z_in_um is the propagation distance


% Ysize = size(ModeFieldsStack, 1) ;
% Xsize = size(ModeFieldsStack, 2) ;
% NumOfModes = size(ModeFieldsStack, 3) ;

    %InitialPhaseStack = angle(ComplexAmplitudeStack) ;
    %A_Stack = abs(ComplexAmplitudeStack) ;
SummedField = sum((ComplexAmplitudeStack.*ModeFieldsStack.*exp(-1i*BetaStack*1e-6*z_in_um)), 3) ;

end
