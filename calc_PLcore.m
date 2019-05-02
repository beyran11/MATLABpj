function PLcore = calc_PLcore(f,Acmag,lc,Bm)

%% Loss method comparison for Metglas POWERLITE C-cores
%% [1]Irma Villar, Unai Viscarret, Ion Etxeberria-Otadui, and Alfred Rufer
% "Global Loss Evaluation Methods for Nonsinusoidally Fed Medium-Frequency 
% Power Transformers", IEEE TRANSACTIONS ON INDUSTRIAL ELECTRONICS,
% VOL. 56, NO. 10, OCTOBER 2009 
%% [2] Irma Villar, Alfred Rufer,Unai Viscarret, Frederic Zurkinden and
% Ion Etxeberria-Otadu, "Analysis of Empirical Core Loss Evaluation
% Methods For Non-Sinusoidally Fed Medium Frequency Power Transformers"
% 2008 IEEE, conf pub


alfa = 1.51 ;
beta = 1.74 ;
K = 6.5 ;
ki = 0.62 ;

Vc = lc*Acmag ;
% Bm = Udc1/(4*f*N1*Acmag) ;

% Calculating Losses density W/kg
Ps_OSE = K*f^alfa*Bm^beta ;
Ps_MSE = (8/pi^2)^(alfa-1)*K*f^alfa*Bm^beta ;
Ps_IGSE = 2^(alfa+beta)*ki*f^alfa*Bm^beta *29.5034e-6 ;

rhoc = 7180 ; % kg/m3

PLcore = Ps_IGSE*rhoc*Vc ;
