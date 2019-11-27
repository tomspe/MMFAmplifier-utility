

function F = CharacteristicFunctionForLP(u, V, l)

%% the characteristic equation is:
%%% for even l :
%%% Jl(u)/Jl+1(u)/u = Kl(w)/Kl+1(w)/w
%%% for odd l :
%%% u*Jl-1(u)/Jl(u) = -w*Kl-1(w)/Kl(w)


ErrorIndices = find( u>=V | u==0 ) ;

w = sqrt(V^2-u.^2) ;
if mod(l,2) == 0
    F = besselj(l,u)./besselj(l+1,u)./u  -  besselk(l,w)./besselk(l+1,w)./w ;
else
    F = u.*besselj(l-1,u)./besselj(l,u)  +  w.*besselk(l-1,w)./besselk(l,w) ;
end
F(ErrorIndices) = Inf ;

%figure ; plot(u, F) ; grid on ; %ylim([-4  4])

end


