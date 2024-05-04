function [lambda_T]=getlambda_ZD(T)
% T is in K

T_c = T - 273.15;
Tg = -113.6;
C1 = 13.62;
C2 = 59.72;
lambda = 5.03e4;
aT0 = 10.^(-C1*(T_c-Tg)./(C2+T_c-Tg));
lambda_T = aT0*lambda;

end