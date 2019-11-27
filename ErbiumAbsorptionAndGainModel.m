

function [ PumpAbsorptionCoefficient, PumpSaturationIntensity, MaximalInversionAchievable, SignalGainCoefficientForTotalInversion, SignalAbsorptionCoefficientFor0Inversion ] = ErbiumAbsorptionAndGainModel(DopingFactor)
h_planck = 6.62e-34 ;
c = 3e8 ;
                                 % for doping factor of 1, the
                                 % concentration is 1000ppm
        ErConcentration = DopingFactor*2.5e25 ; % Ytterbium doping concentration in m^-3 
        

        UpperLevelLifeTime = 10e-3  ;  % in sec  % taken from the table in Barnard IEEE 1994, also from the book of Dutta pg. 57
   %         %****
        % for pumping @ 976nm :
        lambda_p = 980e-9 ;
        CSabsorp_pump = 1.7e-25 ;   %in m^2    %  lowest value in 1991 IEEE paper Barnes Laming et al. (the Ge and AL doped silica)
                                                % see also Dutta pg. 57 and 1994 IEEE Bernard 
        CSemission_pump = 0 ;  %in m^2   % upper level in Erbium is immediately depopulated by SE, so pump is not regenerated as in Yb! 
%         % for pumping @ 915nm :
%         lambda_p = 915e-9 ;
%         CSabsorp_pump = 0.8e-24 ;   %in m^2
%         CSemission_pump = (100/97 - 1)*CSabsorp_pump ;  %in m^2   %because the maximum possible population n2 is 97%   
        %****
        CSemission_sig = 4e-25 ;  %in m^2   % from IEEE paper Barnes Laming et al. (the Ge and AL doped silica)
        CSabsorp_sig = 2.25e-25 ;  %in m^2  

        CSemission_sig = CSemission_sig/2 ;  % hypotheticaly - moving more to the right on the spectrum
        CSabsorp_sig = CSabsorp_sig/2 ;
        
        
%% For the pump:
PumpAbsorptionCoefficient =  CSabsorp_pump*ErConcentration  ; % 1/the typical absorption length of pump 

PumpSaturationIntensity = h_planck*c/lambda_p/(CSabsorp_pump + CSemission_pump)/UpperLevelLifeTime ;  % see point #2 in Paschotta's "useful numbers"

MaximalInversionAchievable = CSabsorp_pump/(CSabsorp_pump + CSemission_pump) ;  % this ratio should be 50% for 976nm and 97% for 915nm


%% For the signal:

SignalGainCoefficient = CSemission_sig*ErConcentration ;
SignalAbsorptionCoefficient = CSabsorp_sig*ErConcentration ;

% see the 2007 IEEE Foster paper
% also matches the 1990 Journal Light. Tech. by Desurvire
 KKfactor = -0.5 ;   
 
  %KKfactor = 2 ;    % hypotheticaly - moving more to the right on the spectrum
     
SignalGainCoefficientForTotalInversion = SignalGainCoefficient + 1i*KKfactor*SignalGainCoefficient ;
SignalAbsorptionCoefficientFor0Inversion = SignalAbsorptionCoefficient + 1i*KKfactor*SignalAbsorptionCoefficient ;
 
 
          % % % based on the lorentzian line simulation:
% RIC_real = 2.2513e-32*YbConcentration ;  % this is the change in n (dimensionless) which represents phase shift
% %
% RIC_real = (7.4/0.22)*RIC_real;   % bending the result to fit the chiang leger paper

end