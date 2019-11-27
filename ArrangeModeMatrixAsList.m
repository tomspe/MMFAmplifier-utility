

function list = ArrangeModeMatrixAsList(Matrix, ValidityMatrix)

% Arrange the values in 'Matrix' (amplitudes/beta values) in a list
% since the entries could be 0 either because the amplitude is 0, or the
% mode does not exist at all, we need the 'ValidityMatrix' - where all the
% existing modes are represented by a non-zero value


% start from the left, work column by column

[ L, M ] = size(Matrix) ;
NumOfSolutions = length(find(ValidityMatrix)) ;
% the first term is for the # of modes in the cos set, the 2nd term is for
% the # of modes in the sin set, which we take to be smaller (we do not
% count the first row twice)
NumOfModes = ( NumOfSolutions ) + ( NumOfSolutions - M ) ;
list = zeros(1, NumOfModes) ;

counter1 = 0 ;
counter2 = NumOfSolutions ;
for c = 1:M
    for r = 1:L
        if ValidityMatrix(r,c) ~= 0
            counter1 = counter1 + 1 ;
            list(counter1) = Matrix(r,c) ;
            if r > 1
                counter2 = counter2 + 1 ;
                list(counter2) = Matrix(r,c) ; 
            end
        end
    end
end

%k = sub2ind([ L, M ], r, c) ;

end
        