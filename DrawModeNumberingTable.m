function DrawModeNumberingTable(ValidityMatrixPump, ValidityMatrixSignal)

%For pump
[ L, M ] = size(ValidityMatrixPump) ;
NumOfSolutions = length(find(ValidityMatrixPump)) ;

Data1 = cell(L, M) ;   Data2 = cell(L, M) ;
counter1 = 1 ;   counter2 = NumOfSolutions + 1 ;
for ColumnIndex = 1:M
    NumOfSolutionsInThisColumn = find(ValidityMatrixPump(:,ColumnIndex), 1, 'last') ;
    for k = 1:NumOfSolutionsInThisColumn
        Data1{k,ColumnIndex} = counter1 ;
        counter1 = counter1 + 1 ;
        if k > 1
            Data2{k,ColumnIndex} = counter2 ;
            counter2 = counter2 + 1 ;
        end
    end
end
        

%Same thing - For signal
[ L, M ] = size(ValidityMatrixSignal) ;
NumOfSolutions = length(find(ValidityMatrixSignal)) ;

Data3 = cell(L, M) ;   Data4 = cell(L, M) ;
counter1 = 1 ;   counter2 = NumOfSolutions + 1 ;
for ColumnIndex = 1:M
    NumOfSolutionsInThisColumn = find(ValidityMatrixSignal(:,ColumnIndex), 1, 'last') ;
    for k = 1:NumOfSolutionsInThisColumn
        Data3{k,ColumnIndex} = counter1 ;
        counter1 = counter1 + 1 ;
        if k > 1
            Data4{k,ColumnIndex} = counter2 ;
            counter2 = counter2 + 1 ;
        end
    end
end
        
ModeNumbering(Data1, Data2, Data3, Data4) ;

end

