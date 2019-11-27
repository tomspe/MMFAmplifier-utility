
function u_values = Find_u_Solutions_by_Mitschke(V, tolerance_epsilon)

% The function finds LP modes by solving for u, in the method outlined by
% Fedor Mitschke's textbook: in pairs of the l index
% for example all LP0x and LP1x modes are alternating roots 
% in the same w vs. u graph, then all roots of LP2x and LP3x similarly are
% found as alternating roots in the next order, etc.


u_values = zeros(100) ;  % we arrange this solution matrix thus:
% the 1st indice (the row #) is l, the angular dependence
% the 2nd indice (the column #) is m, the radial dependence
% we go over l from 0 to l_max where Jl no longer has zeros in the 0-V range, 
% and for each such l find all m values - by counting how many alternating 
% ranges between roots are found in that range


l = 0 ;    % azimuthal index for going over bessel orders (0, 2, 4, ..)
NumOfZerosCeil = round(V/2)  ;   % initial guess of # of zeros in the 0-V range exaggerated upwards
GoOn = 1 ;
while GoOn
    % first we find all bessel zeros, this defines the intervals
    EvenZeros = besselzero(l, NumOfZerosCeil, 1) ;
    OddZeros = besselzero(l+1, NumOfZerosCeil, 1) ;
    
    % the guessed number of solutions should be larger than the true number
    % we find the first root of the even Bessel func that's beyond V
    StopInd = find(EvenZeros > V, 1, 'first') ;
    if StopInd == 1
        GoOn = 0 ;
    else
        for m = 1:(StopInd-1)
            % find solution LPlm in the interval below the zero of the even bessel func
            if m == 1    % special case of leftmost interval, left border is simply u=0
                LeftLimit = 0 ;
            else
                LeftLimit = OddZeros(m-1) ;
            end
            u0 = (LeftLimit + min(V, EvenZeros(m)))/2 ;    % the middle of the range - starting point for search of zero
            u_solution = fzero(@(u)  CharacteristicFunctionForLP(u, V, l), u0) ;
            if CharacteristicFunctionForLP(u_solution, V, l) <= tolerance_epsilon   % verify - in case the fzero function did not converge, do nothing
                u_values(l+1,m) = u_solution ;
            end
                
            % find solution LPlm in the interval above the zero of the even bessel func
            u0 = (EvenZeros(m) + min(V, OddZeros(m)))/2 ;    % the middle of the range - starting point for search of zero
            u_solution = fzero(@(u)  CharacteristicFunctionForLP(u, V, l+1), u0) ;
            if CharacteristicFunctionForLP(u_solution, V, l) <= tolerance_epsilon    % verify - in case the fzero function did not converge, do nothing
                u_values(l+2,m) = u_solution ;
            end 
            disp([ 'u values for LP' num2str(l) num2str(m) ', LP' num2str(l+1) num2str(m) ' - done' ]) 
        end
    end
    
    l = l + 2 ;   % jump up 2 because each iteration solves 2 intervals, the even and the odd
end

l_max = find(u_values(:,1), 1, 'last') ;
m_max = find(u_values(1,:), 1, 'last') ;
u_values = u_values(1:l_max, 1:m_max) ;  % crop matrix by highest index reached


end