
% Correction - Tij and Tji should not be conjugates, but identical
function T = ArrangeListAsSymmetricMatrix(vec, N)

% vec should be arranged into the upper triangle + the diagonal of a N*N
% matrice, where N is the number of modes, therefore it should have length
% of (N+1)*N/2. bottom triangle is then complemented by Hermitian
% characteristic: Tij = Tji* for all i not equal to j (hence the "almost" Hermitian)

if length(vec) ~= (N + 1)*N/2
    errdlg('Error in # of modes!') ;
    T = [] ;
    return
end


Mat = zeros(N) ;
for r = 1:N
    ei = r*N - (r - 1)*r/2 ;
    si = ei - (N - r) ;
    Mat(r, r:N) = vec(si:ei) ; 
end
% Mat is now a matrix with a non-zero diagonal and upper triangle
% we want to duplicate the upper triangle to the lower with conjugation,
% without touching the diagonal

Mat_wo_diag = Mat - diag(diag(Mat)) ;
T = Mat_wo_diag + Mat_wo_diag.' + diag(diag(Mat)) ;

end