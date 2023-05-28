function Xdot=odefcn_endogAlpha(t,X,year,Epath,F_nonco2,J,G,CS)
%CS climatesensitivity;
%J,G are 1 dimentional arrays with parameters from Joos and Geoffroy respectively
%X contains 4 carbon boxes in GtCO2 above preindustrial (slowest decaying
%first; Temp_atmosphere; Temp deep ocean; cumulative emissions since
%preindustrial
Epath=interp1(year,Epath,t);
F_nonco2=interp1(year,F_nonco2,t);
Alpha=0.0106774818735627*exp(0.0866328405015618*(34.4+0.019*(X(7)-X(1)-X(2)-X(3)-X(4))/3.663+4.165*(X(5)))); 
%coefficients from xfit_exp in SolveAlpha, from FAIR paper 
%571 cumulative emissions in 2015 in GtC
%subtract (X(1)+X(2) X(3) X(4))/3.666 from cumulative emissions according to paper FAIR
Xdot= [J(1)*Epath; ...
       J(2)*Epath-J(5)/Alpha*X(2); ...
       J(3)*Epath-J(6)/Alpha*X(3); ...
       J(4)*Epath-J(7)/Alpha*X(4); ...
       (CS*G(3)*log(1+(X(1)+X(2)+X(3)+X(4))/(3.663*588))/log(2)+F_nonco2)/G(1) ... % 588 preindustrial atmosrpheric concentration in GtC (DICE)
       - G(3)/G(1)*X(5)- G(4)/G(1)*(X(5)-X(6));...
       G(4)/G(2)*(X(5)-X(6));...
       Epath];
end