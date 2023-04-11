% 0.2 HP 1500 rpm, T=2.481 N.m, 1.3 amp

Vrms=120/sqrt(3); 
Vmax=sqrt(2)*Vrms;  
fe  =60; we=2*pi*fe;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          % system frequency
J   =0.0068;              
P   = 4;                 

ras=14; rbs=14; rcs=14;             
rar=7.7; rbr=7.7; rcr=7.7;
Lls=9/we; Lms=155/we; Llr=9/we;           

Tm=2.418;

V1s=Vmax;
V3s=0*Vmax/3;
V5s=0*Vmax/5;
V7s=0*Vmax/7;
V11s=0*Vmax/11;
V13s=0*Vmax/13;
V17s=0*Vmax/17;
V19s=0*Vmax/19;
V23s=0*Vmax/23;
V25s=0*Vmax/25;
V29s=0*Vmax/29;
% 
V1r=0*Vmax/3;
V5r=0*V1r/5;
V7r=0*V1r/7;
V11r=0*V1r/11;
V13r=0*V1r/13;
V17r=0*V1r/17;
V19r=0*V1r/19;
V23r=0*V1r/23;
V25r=0*V1r/25;
V29r=0*V1r/29;

%%%%%%%%%%%%%%%% TODO CORRECTO %%%%%%%%%%%%%%%%%%%%

% Vrms=220/sqrt(3); Vr=50;          
% Vmax=sqrt(2)*Vrms;
% Vrh=Vmax/Vr;
% fe  =50; we=2*pi*fe;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          % system frequency
% J   =0.0068;              % inertia kg.m2 
% P   = 4;                  % number of poles
% 
% ras=45;%14;               % stator resistences ohms
% rbs=45;%14;
% rcs=45;%14;
% rar=38;%7.7;              % rotor resistences ohms
% rbr=38;%7.7;
% rcr=38;%7.7;
% Lls=38;%9/we;             % stator inductance
% Lms=700;%155/we;          % inductance
% Llr=38;%9/we;             % rotor inductance
% 
% V1s=Vmax;
% V3s=0*Vmax/3;
% V5s=0*Vmax/5;
% V7s=0*Vmax/7;
% 
% V1r=0*Vmax/Vrh;
% V3r=0*V1r/3;
% V5r=0*V1r/5;
% V7r=0*V1r/7;
% 
% Tm=0;%2.418;
