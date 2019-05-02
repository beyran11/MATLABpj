function Ls = calc_Ls(N1,N2,Np,A,B,C,D,M,df,do1,do2,di1,di2,dg,mu0)
% leakage inductance is given for the concentric (16)
% or Horizontal
concentric = 0;
if concentric
    
    hwmax = C-2*df ;
    
    m1 = (N1+1)*do1/hwmax ;
    m2 = (N2+1)*do2/hwmax ;
    hw = (N1+1)/m1*do1 ;
    dw = m1*do1 + m2*do2 + (m1-1)*di1 + (m2-1)*di2 + M*dg ;
    lw = 2*(2*A+Np*D)+2*pi*(dw/2+df) ;
    Ls = mu0*(N1/M)^2*lw/hw*(m1*do1/3+m2*do2/3+M*dg) ;
end

% for the split windings (17)
split = 1;
if split
    dair_main = 0 ;
    hwmax = B-dair_main-df ;
    
    m1 = N1*(do1+di1)/(hwmax-di1) ;
    m2 = N2*(do2+di2)/(hwmax-di2) ;
    hw = (N1/m1)*do1+(N1/m1-1)*di1 ;
    dw = (m1+1)*do1 + (m2+1)*do2 + M*dg ;
    lw = 2*(2*A+Np*D)+2*pi*(dw/2+df) ;
    
    Ls = mu0*(N1/M)^2*lw/hw*(m1*do1/3+m2*do2/3+...
        di1*(m1-1)/(2*m1)+di2*(m2-1)/(2*m2)+M*dg) ;
end