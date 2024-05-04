function [lambda_T]=getlambda_ZT(T)
% T is in K

T_c = T - 273.15;
Tg = -112.2;
C1 = 23.22;
C2 = 45.81;
lambda = 4.02e13;
aT0 = 10.^(-C1*(T_c-Tg)./(C2+T_c-Tg));
lambda_T = aT0*lambda;

end