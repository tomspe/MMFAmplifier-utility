

function [ PumpAbsorptionCoefficient, PumpSaturationIntensity, MaximalInversionAchievable, SignalGainCoefficientForTotalInversion, SignalAbsorptionCoefficientFor0Inversion ] = YtterbiumAbsorptionAndGainModel(DopingFactor)
h_planck = 6.62e-34 ;
c = 3e8 ;
                                 % for doping factor of 1, the
                                 % concentration is 1000ppm
        YbConcentration = DopingFactor*1e25 ; % Ytterbium doping concentration in m^-3 
        
        lambda_s = 1030e-9 ;
        UpperLevelLifeTime = 0.8e-3  ;  % in sec
        %****
        % for pumping @ 976nm :
        lambda_p = 976e-9 ;
        CSabsorp_pump = 2.6e-24 ;   %in m^2
        CSemission_pump = CSabsorp_pump ;  %in m^2 
%         % for pumping @ 915nm :
%         lambda_p = 915e-9 ;
%         CSabsorp_pump = 0.8e-24 ;   %in m^2
%         CSemission_pump = (100/97 - 1)*CSabsorp_pump ;  %in m^2   %because the maximum possible population n2 is 97%   
        %****
        CSemission_sig = 6.2734e-25 ;  %in m^2
        CSabsorp_sig = 45e-27 ;  %in m^2  % from the 1997 paschotta paper fig. 2

%% For the pump:
PumpAbsorptionCoefficient =  CSabsorp_pump*YbConcentration  ; % 1/the typical absorption length of pump 

PumpSaturationIntensity = h_planck*c/lambda_p/(CSabsorp_pump + CSemission_pump)/UpperLevelLifeTime ;  % see point #2 in Paschotta's "useful numbers"

MaximalInversionAchievable = CSabsorp_pump/(CSabsorp_pump + CSemission_pump) ;  % this ratio should be 50% for 976nm and 97% for 915nm


%% For the signal:
% RIC_imag_gain = (1/2)*CSemission_sig*YbConcentration*(lambda_s/2/pi) ;  % this is the change in n (dimensionless) which represents gain for the signal
% RIC_imag_absorp = -(1/2)*CSabsorp_sig*YbConcentration*(lambda_s/2/pi) ;
%               % factor 1/2 in the exponent of the field, so that alpha = sigma*N, where alpha is the coefficient for intensity 
%               
% % see the chiang leger paper
% KKfactor = 7.4 ;
% SignalComplexRefractiveIndexForTotalInversion = (RIC_real_gain + 1i*KKfactor*RIC_imag_gain) ;
% SignalComplexRefractiveIndexFor0Inversion = (RIC_real_absorp + 1i*KKfactor*RIC_imag_absorp) ;

SignalGainCoefficient = CSemission_sig*YbConcentration ;
SignalAbsorptionCoefficient = CSabsorp_sig*YbConcentration ;
% see the chiang leger paper
 KKfactor = 7.4 ;   
SignalGainCoefficientForTotalInversion = SignalGainCoefficient + 1i*KKfactor*SignalGainCoefficient ;
SignalAbsorptionCoefficientFor0Inversion = SignalAbsorptionCoefficient + 1i*KKfactor*SignalAbsorptionCoefficient ;
 
 
          % % % based on the lorentzian line simulation:
% RIC_real = 2.2513e-32*YbConcentration ;  % this is the change in n (dimensionless) which represents phase shift
% %
% RIC_real = (7.4/0.22)*RIC_real;   % bending the result to fit the chiang leger paper

end