function Tdot=odefcn_inertia(t,T,Forcing,year,G)
Forcing=interp1(year,Forcing,t);
Tdot=[(Forcing - G(3)*T(1)- G(4)*(T(1)-T(2)))/G(1); ...
    G(4)*(T(1)-T(2))/G(2)];
end