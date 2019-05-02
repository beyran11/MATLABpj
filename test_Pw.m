% function calc_PLw
clc
clear all
close all


m = 5 ;
h = 1 ;
f = 2e3;
mu0 = 4*pi*1e-7 ;
s_cu = 5.688e7 ;
ds = 0.2e-3 ; % 0.1<ds<0.355
% dw = sqrt(pi/4)*ds
ns = 600 ; % 3<ds<2520
depth = 1/sqrt(pi*h*f*mu0*s_cu) 
delta = ds/depth 
Zeta1 = (sinh(2*delta)+sin(2*delta))/(cosh(2*delta)-cos(2*delta)) ;
Zeta2 = (sinh(delta)-sin(delta))/(cosh(delta)+cos(delta)) ;
Fr = delta*(Zeta1+2/3*(m^2-1)*Zeta2);

omga1 = 0;
omga2 = 0;
gama = 0.1309 ;
Udc1 = 1100 ;
Udc2 = 1100 ;
Ls = 0.480e-3;
Uac1 = Udc1*2./(h*pi).*cos(h*omga1).*(1-cos(h*pi));
Uac2 = Udc2*2./(h*pi).*cos(h*omga2).*(1-cos(h*pi));
dU = sqrt(Uac1.^2 + Uac2.^2 - 2*Uac1.*Uac2.*cos(h*gama));
Iac = dU./(2*pi*f*h*Ls) ;

lw = 684*1e-3 ;
twisting_factor = 1.05 ;
% N1*lw*twisting_factor ;
RDC = 27.4e-3 ;

Rac = RDC*Fr;
PLac = RDC*Fr*Iac^2 ;

%% Round Litz wires skin effect losses


gama = ds/(depth*sqrt(2)) ;
ls = 1 ;
% Rdcs = ls*4/(pi*s_cu*ds^2) ;
% [Lammeraner and Štafl, 1966]
taw_1 = 2/gama + gama^3/96 
% Pskin_Lam = Rdcs*gama/4*taw_1*3.8 ;
% [Lotfi and Lee, 1993]
% Wave number
alfa = (1+1j)/depth;
I0 = besseli(0,alfa*ds/2);
I1 = besseli(1,alfa*ds/2);
taw_1 = 1/sqrt(2)*real((1+1j)*I0/I1);

% Pskin_Lot = Rdcs*gama/4*taw_1*3.8 ; 

%% Round Litz wires Internal Proximity Effect Losses 
db = (ds*1.1)*sqrt(2/pi*sqrt(3)*ns)
% [Lammeraner and Štafl, 1966]
taw_2 = -gama^3/16
PLinternal = -(2*pi*gama/s_cu)*ns/(8*pi^2*(db/2)^2)*taw_2*120
% [Lotfi and Lee, 1993]
% taw_2 = 1/sqrt(2)*real((1-1j)*I0*conj(I1))/abs(I0)^2
% PLinternal_Lot = -(2*pi*gama/s_cu)*ns/(8*pi^2*(db/2)^2)*taw_2*120

%% Round Litz wires External Proximity Effect Losses 
He = 1e3;
PLexternal = -(2*pi*gama/s_cu)*ns*He^2*taw_2

%% total loss
%  [Bartoli et al., 1996]
tb = 1.1*db ;
ts = 1.1*ds ;
eta1 = ds/tb*sqrt(pi/4);
eta2 = ds/ts*sqrt(pi/4);
pf = ns*(ds/db)^2 ;
PLw = Iac^2/2*RDC*gama/2*(1/ns*taw_1 - ...
    2*pi*(4*(m^2-1)/3+1)*ns*(eta1^2+eta2^2*pf/(2*pi*ns))*taw_2)

%  [Lammeraner and Štafl, 1966]
PLw = Iac^2*RDC*(1+gama^4/192*(1/6+pi^2*ns*pf/4*(16*m^2-1+24/pi^2)))
