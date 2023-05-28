function Mdot=odefcn_decay(t,M,Epath,year,J,Alpha)
Epath=interp1(year,Epath,t);
Mdot= [J(1)*Epath;  J(2)*Epath-J(5)/Alpha*M(2); ...
                    J(3)*Epath-J(6)/Alpha*M(3); ...
                    J(4)*Epath-J(7)/Alpha*M(4)];

end
