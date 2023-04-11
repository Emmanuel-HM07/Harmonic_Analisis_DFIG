function [A,B,U]=parameter_induction_machine(t,X)

% Reference: P.C. Krause, et. al, Analysis of Electrical Machinery and
% Drive Systems, second Edition, IEEE Press, Wiley-Interscience, Chapter 4.
% dIsr=AIrs+BVsr induction machine state equation with rotor variables refered to
% to stator
% Isr=[ias ibs ics iar ibr icr]'
% Vsr=[vas vbs vcs var vbr vcr]'
% dwr=(P/2J)(Te-Tm)
% dOr=wr;
% X=[Isr wr Or]
%

23hp_induction_machine;

% 3 HP induction machine paramenters:

% state variables

ias=X(1); ibs=X(2); ics=X(3);
iar=X(4); ibr=X(5); icr=X(6);
wr =X(7); Or =X(8);

Iabcs=[ias;ibs;ics];
Iabcr=[iar;ibr;icr];
Isr=[Iabcs;Iabcr];

% stator voltage and rotor voltage

Vmax=sqrt(2)*Vrms;

% vas=Vmax*cos(we*t);
% vbs=Vmax*cos(we*t-2*pi/3);
% vcs=Vmax*cos(we*t+2*pi/3);

vas=Vmax*sin(we*t);
vbs=Vmax*sin(we*t-2*pi/3);
vcs=Vmax*sin(we*t+2*pi/3);

var=0;
vbr=0;
vcr=0;

vabcs=[vas;vbs;vcs];
vabcr=[var;vbr;vcr];

Vsr=[vabcs;vabcr];

% Induction machine matrices of parameters

Rs=[ras 0   0
    0   rbs 0
    0   0   rcs];

Rr=[rar 0   0
    0   rbr 0
    0   0   rcr];

Ls=[ Lls+Lms -1/2*Lms -1/2*Lms
    -1/2*Lms  Lls+Lms -1/2*Lms
    -1/2*Lms -1/2*Lms  Lls+Lms ];

Lr=[ Llr+Lms -1/2*Lms -1/2*Lms
    -1/2*Lms  Llr+Lms -1/2*Lms
    -1/2*Lms -1/2*Lms  Llr+Lms ];

Lsr=Lms*[cos(Or)        cos(Or+2*pi/3) cos(Or-2*pi/3)
         cos(Or-2*pi/3) cos(Or)        cos(Or+2*pi/3)
         cos(Or+2*pi/3) cos(Or-2*pi/3) cos(Or)       ];
     
pLsr=-wr*Lms*[sin(Or)        sin(Or+2*pi/3) sin(Or-2*pi/3)
              sin(Or-2*pi/3) sin(Or)        sin(Or+2*pi/3)
              sin(Or+2*pi/3) sin(Or-2*pi/3) sin(Or)       ];
L=[Ls   Lsr
   Lsr' Lr];

R=[Rs    pLsr
   pLsr' Rr];
          
B=inv(L);
A=-B*R;

dIsr=A*Isr+B*Vsr;

% Torque and dwr

i1=ias*(iar-ibr/2-icr/2)+ibs*(ibr-iar/2-icr/2)+ics*(icr-ibr/2-iar/2);
i2=ias*(ibr-icr)+ibs*(icr-iar)+ics*(iar-ibr);
Te=-(P/2)*Lms*(i1*sin(Or)+sqrt(3)/2*i2*cos(Or));

dwr=P*(Te-Tm)/(2*J);
dOr=wr;

% state variables

dX=[dIsr;dwr;dOr];
          
end