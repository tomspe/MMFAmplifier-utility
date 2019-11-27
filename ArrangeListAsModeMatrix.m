

function Matrix = ArrangeListAsModeMatrix(List, ValidityMatrix, FirstRowIsNullFlag)

% Arrange the values in the list back to matrix form
% since the entries could be 0 either because the amplitude is 0, or the
% mode does not exist at all, we need the 'ValidityMatrix' - where all the
% existing modes are represented by a non-zero value

% the FirstRowIsNullFlag is 'on' when we are dealing with the 2nd set of
% degenerate modes and the list is shorter (does not contain values for the
% 1st row of thz table, where modes are not degenerate)

% start from the left, work column by column 
[ L, M ] = size(ValidityMatrix) ;
cutoff_indices = find(ValidityMatrix == 0) ;
NumOfModes = L*M - length(cutoff_indices) ;
% check that # of modes is the same as in the list
if (~FirstRowIsNullFlag && NumOfModes ~= length(List)) || (FirstRowIsNullFlag && (NumOfModes - M) ~= length(List)) 
    errordlg('No agreement in the Number of Modes!', 'Modal') ;
    Matrix = [] ;
    return
end

Matrix = zeros(L, M) ;
counter = 1 ;
for k = 1:L*M
    if ValidityMatrix(k) ~= 0 
        [r,~] = ind2sub([ L, M ], k) ;
        if FirstRowIsNullFlag && r == 1 
            % fill in the first row with 0 since it is null anyway
            Matrix(k) = 0 ;
        else
            Matrix(k) = List(counter) ;
            counter = counter + 1 ;
        end
    end
end


end
        