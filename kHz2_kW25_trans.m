clc
clear all
close all


m = 5 ;
f = 2e3 ;
mu0 = 4*pi*1e-7 ;

ds1 = 0.2e-3 ; % 0.1<ds<0.355
ds2 = ds1 ;
ns1 = 600 ; % 3<ds<2520
ns2 = ns1 ;
gama = 0.1309 ;
Udc1 = 1100 ;
Ls = 0.480e-3 ;

Bmag = 0.5 ;
lgap = 0.5e-3 ;

%% Control Variables
Np = 7 ;
M = 1 ;
% Bm = 0.23 ;
dg = 1.5e-3 ;


m1 = 5 ;
m2 = 5 ;

%% secifications

db1 = 1.1*ds1*sqrt(2*sqrt(3)*ns1/pi) ;
db2 = 1.1*ds2*sqrt(2*sqrt(3)*ns2/pi) ;

pf1 = ns1*(ds1/db1)^2 ;
pf2 = ns2*(ds2/db2)^2 ;
%% constants

mur = 1000 ;
rhoc = 7180 ;
d_former = 6e-3 ;
d_intra = 1e-3 ;
Lf = 0.84 ;  % lamination factor lf is 0.84

%% Split windings
NC = 2 ;
d_iso = 5e-3;
d_air = 7e-3;
d_former = 6e-3;
d_intra = 1e-3;
rho_copper = 8960 ; % density of copper is 8.96 g/cm3

%% Shell-Type Split Winding

%% Villar Thesis 2010
A = 25.8*1e-3 ;
% B = 67.0*1e-3 ;
% C = 97.8*1e-3 ;
D = 25.0*1e-3 ;

Ac = 2*Np*A*D ;
Acmag = Ac*Lf ;
N1 = Udc1/(4*f*Bmag*Acmag) 
% N1 = 38 ;
N2 = N1 ;

w_w1 = m1*db1 ;
w_w2 = m2*db2 ;
w_w = w_w1+w_w2+d_iso ;
h_w1 = (N1/m1+1)*db1 + (N1/m1-1)*d_intra ;
h_w2 = (N2/m2+1)*db1 + (N2/m2-1)*d_intra ;
l_w1 = (d_former + w_w1)*2*pi + 2*(Np*D+2*A) ;
l_w2 = (d_former + w_w2)*2*pi + 2*(Np*D+2*A) ;

C = 2*d_former + w_w1 + w_w2 + d_iso ;
B = d_former + max(h_w1,h_w2) + d_air ;
E = B + 2*A ;
F = C + 2*A ;
lc = 2*(B+C)+4*A 
lmag = lc ; 

%% Leakage inductance Ls calculation
Ls = calc_Ls(N1,N2,Np,A,B,C,D,M,d_former,db1,db2,d_intra,d_intra,dg,mu0)

%% Magnetization inductance Lm calculation
Lm = calc_Lm(Acmag,N1,mu0,lgap,lmag,mur) 
%% Core Losses calculation
% PLcore = calc_PLcore(Udc1,f,N1,Acmag,lc);
PLcore = calc_PLcore(f,Acmag,lc,Bmag)
%% Winding Losses calculation
PLw_p =calc_PLw(ds1,ns1,m1,f,N1,l_w1) ;
PLw_s =calc_PLw(ds2,ns2,m2,f,N2,l_w2) ;
Ploss = PLcore + PLw_p + PLw_s

%% Copper weight
Vw = l_w1*h_w1*w_w1 + l_w2*h_w2*w_w2 ;
Weight_w = Vw*rho_copper ;

%% Metglas Wolume
Vc = lc*Ac ;
Wc = rhoc*Vc ;

%% Final Volume
Vt =2*E*F*(Np*D + d_former + max(h_w1,h_w2)) ;

