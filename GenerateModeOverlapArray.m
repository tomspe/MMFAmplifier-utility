
function AA = GenerateModeOverlapArray(u_matrix, ModeStack)


% for N supported modes, we prepare (N+1)*N/2 overlap matrices,
% one for each element in the upper triangle (+ the diagonal) of the NxN matrix
% describing all possible mode couplings (it's Hermitian so there is no need to calculate the bottom triangle) 

NumOfSpatialPoints = size(ModeStack, 1) ;
N = size(ModeStack, 3) ;  
AA = zeros(NumOfSpatialPoints, NumOfSpatialPoints, N*(N + 1)/2) ;
c = 0 ;
for i = 1:N
    for j = i:N
        c = c + 1 ;
        AA(:,:, c) = ModeStack(:,:,i).*ModeStack(:,:,j) ;
    end
end

if c ~= N*(N + 1)/2
    errordlg('Error in creation of mode couplings array!')
end

end