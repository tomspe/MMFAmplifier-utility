

function ModeFieldsStack = GenerateStackOfModeFields(u_matrix, V, XaxisExtended, CropMar)

debug = 0 ;
% the 0 entries in the u value matrix are non-existent modes (beyond
% cut-off) we count how many modes exist
[ L, M ] = size(u_matrix) ;
cutoff_indices = find(u_matrix == 0) ;
NumOfModes = L*M - length(cutoff_indices) ;

w_matrix = sqrt(V^2 - u_matrix.^2) ;

% we arrange the stack by going over the modes in columns,
% so first all the LPx1 modes are computed, then all LPx2, etc.
% *note - this is contrary to the order in which we find them by solving
% the characteristic equation, the latter process is performed by rows
StackExtended = zeros(length(XaxisExtended), length( XaxisExtended),2*NumOfModes) ;   % factor 2 for theta degeneracy
                                                                    % "Dead modes" (for those that are not degenerat) will be cropped out later to not waste space        

% preparing the spatial coordinates - the axes were prepared so that the
% core radius is 1
[ XX, YY ] = meshgrid(XaxisExtended, XaxisExtended) ; 
r = sqrt(XX.^2 + YY.^2) ; 
theta = atan2(YY, XX) ;
Core = r<= 1 ;
Clad = r > 1 ;
c = 0 ;         %counter for set1
CounterForModeSet2 = 0 ;   % counter of modes in set 2 (where some are 'dead'
for k = 1:L*M
    if u_matrix(k) ~= 0
        c = c + 1 ;
        [ row, col ] = ind2sub([ L M ], k) ;
        l = row - 1 ;      % the angular index l starts from 0 not from 1, as opposed to the radial index m (LP01 is the first mode!)
        m = col ;   % not actually used in the formula for the fields
        
        OnlyRadialComp_core = Core.*besselj(l, u_matrix(k)*r) ;
        OnlyRadialComp_clad = Clad.*besselk(l, w_matrix(k)*r)*(besselj(l,u_matrix(k))/besselk(l,w_matrix(k))) ;
        Fcore1 = OnlyRadialComp_core.*cos(l*theta) ;
        Fcore2 = OnlyRadialComp_core.*sin(l*theta) ;
        Fclad1 = OnlyRadialComp_clad.*cos(l*theta) ;
        Fclad2 = OnlyRadialComp_clad.*sin(l*theta) ;
        % for normalization:
        Pcore1 = sum(sum(Fcore1.^2)) ; %/sum(sum(Core)) ;
        Pclad1 = sum(sum(Fclad1.^2)) ; %/sum(sum(Clad)) ;
            Pcore2 = sum(sum(Fcore2.^2)) ; %/sum(sum(Core)) ;
            Pclad2 = sum(sum(Fclad2.^2)) ; %/sum(sum(Clad)) ;
        % saving field to the stack, normalized to power=1
        StackExtended(:,:,c) = (Fcore1 + Fclad1)/sqrt(Pcore1 + Pclad1) ;
        if l > 0
            CounterForModeSet2 = CounterForModeSet2 + 1 ;
            StackExtended(:,:,NumOfModes + CounterForModeSet2) = (Fcore2 + Fclad2)/sqrt(Pcore2 + Pclad2) ;
        end

        if debug
            F = Fcore + Fclad ;
            %% power on log scale
           % figure ; imagesc(log(1+ F.^2)) ; colorbar ; 
           % hold on ; plot(0.5+Xsize, 0.5+1.5*Ysize, '*k') ; plot(0.5+2*Xsize, 0.5+1.5*Ysize, '*k') ;  plot(0.5+1.5*Xsize, 0.5+Ysize, '*k') ;  plot(0.5+1.5*Xsize, 0.5+2*Ysize, '*k') ;    
           
%             %% derivatives:
%            figure ;
%            subplot(1,2,1) ; imagesc(F(:,2:end) - F(:,1:(end-1))) ; colorbar ;
%            hold on ; plot(0.5+Xsize, 0.5+1.5*Ysize, '*k') ; plot(0.5+2*Xsize, 0.5+1.5*Ysize, '*k') ;  plot(0.5+1.5*Xsize, 0.5+Ysize, '*k') ;  plot(0.5+1.5*Xsize, 0.5+2*Ysize, '*k') ;
%            subplot(1,2,2) ; imagesc(F(2:end,:) - F(1:(end-1),:)) ; colorbar ; 
%            hold on ; plot(0.5+Xsize, 0.5+1.5*Ysize, '*k') ; plot(0.5+2*Xsize, 0.5+1.5*Ysize, '*k') ;  plot(0.5+1.5*Xsize, 0.5+Ysize, '*k') ;  plot(0.5+1.5*Xsize, 0.5+2*Ysize, '*k') ;
          
           confinement(c) = 100*Pcore/(Pcore + Pclad) ;
           temp = ((u_matrix(k)/V)^2)*(1 - (besselk(l, w_matrix(k)))^2/besselk(l+1, w_matrix(k))/besselk(l-1, w_matrix(k))) ;
           confinement_theory(c) = 100*(1 - temp) ;
        end
        
    end
end

% we crop the margins far from the fiber core
% Also - beacuse of the 'dead' modes, we have some unused layers at the bottom of
% our 3D array, we crop them
ModeFieldsStack = StackExtended((CropMar+1):(end-CropMar),(CropMar+1):(end-CropMar),1:(NumOfModes + CounterForModeSet2)) ;


if debug
    figure ; plot(confinement) ; hold on ; plot(confinement_theory, 'r') ; grid on ; legend('Numerically Integrated', 'Theoretical formula')
end

% checking orthogonality:
%sum(sum(Core.*ModeFieldsStack(:,:,1).*conj(ModeFieldsStack(:,:,2))))/sum(sum(Core))
% Xaxis = XX(1,:) ;
% Yaxis = YY(:,1) ;
 
end