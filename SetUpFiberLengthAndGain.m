
% small auxilary function for retrieving the fiber's Z axis, doping, and
% normalizing the pump amplitudes to the real power in W - before calling any calculation function
function [  Zaxis, DopingFactor, PumpPower ] = SetUpFiberLengthAndGain(handles) 

FiberLength = 1e4*str2double(get(handles.FiberLength_et, 'String')) ;   % convert to um
PumpPower = 1e-3*str2double(get(handles.PumpPower_et, 'String')) ;   % convert to W
DopingFactor = str2double(get(handles.GainFactor_et, 'String')) ;

% allocating the Z axis: L is the total length, dz1 is the step of the propagator - corrected for the fact the original dz is not a perfect divisor of L 
Nsteps = ceil(FiberLength/handles.dz) ; dz1 = FiberLength/Nsteps ;  Zaxis = 0:dz1:FiberLength ;    % in um

end