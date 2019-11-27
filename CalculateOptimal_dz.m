
function dephasing_length = CalculateOptimal_dz(beta_list, amplitude_list, DephasingSamplingFactorAlongZ) 

% the lists correspond to the modes, we'd like to take into account only
% the betas of the modes that have non-zero amplitude, since the less modes we consider - 
% the larger the dephasing length will be 

Betas = beta_list( amplitude_list ~= 0   ) ;

% if length(Betas) > 1
    MaximalBetaDelta = max(Betas) - min(Betas) ;
    dephasing_length = pi/DephasingSamplingFactorAlongZ/MaximalBetaDelta ;
% else
%     dephasing_length = 50e-6 ;
% end

dephasing_length = min(100e-6, dephasing_length) 

%1e6*dephasing_length

end