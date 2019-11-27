

function Stack = DuplicateValuesInListToStack(list, Ysize, Xsize)

% input is a vector where each single entry is a value for a mode
% (Amplitude/Beta), the function returns a XxYxNumOfModes array where each
% XxY layer is a uniform matrix, all entries are the corresponding value


% make sure list is a row not a column
if size(list,1) ~= 1
    list = list.' ;
end
NumOfModes = length(list) ;


%duplicate each value to a column of length XxY 
M = repmat(list, Ysize*Xsize, 1) ;
Stack = reshape(M, Ysize, Xsize, NumOfModes) ;

end