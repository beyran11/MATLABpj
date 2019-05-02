function Lm = calc_Lm(Ac,N1,mu0,lgap,lmag,mur)

% Equation (18) of garcia paper
% Lm = mu0*N1^2*Ac*1e-3/(2*lgap+lmag/mur) ;

Lm = mu0*N1^2*Ac*1e3/(lgap+lmag/mur) ;