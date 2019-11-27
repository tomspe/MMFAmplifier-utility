


function beta_matrix = CalculateBetaFor_u_Solutions(u_matrix, n_core, a, lambda)
% beta is the propagation constant for each mode
%  ~exp(-j*beta*z)

% dephasing length is an aaproximate measure of the z distance over which
% the 2 most different betas acquire a significant phase gap, i.e. the maximal distance that can be
% used in the simulation as an incerement in the Z axis 

% the 0 entries in the u value matrix are non-existent modes (beyond
% cut-off) we count how many modes exist
[ L, M ] = size(u_matrix) ;
cutoff_indices = find(u_matrix == 0) ;
existing_indices = find(u_matrix ~= 0) ;
    %NumOfModes = L*M - length(cutoff_indices) ;


beta_matrix = sqrt(  (n_core*2*pi/lambda)^2 - (u_matrix/a).^2  ) ;
beta_matrix(cutoff_indices) = 0 ;
    % avg_beta = sum(sum(beta_matrix))/NumOfModes ;


end