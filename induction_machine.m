function dX=induction_machine(t,X)

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

induction_machine_23HP; %parameters of the induction machine 

% state variables

ias=X(1); ibs=X(2); ics=X(3);
iar=X(4); ibr=X(5); icr=X(6);
wr =X(7); Or =X(8);

Iabcs=[ias;ibs;ics];
Iabcr=[iar;ibr;icr];
Isr=[Iabcs;Iabcr];

% vas=V1s*cos(we*t)       +V7s*cos(7*we*t)       +V5s*cos(5*we*t)       +V3s*cos(3*we*t); 
% vbs=V1s*cos(we*t-2*pi/3)+V7s*cos(7*we*t-2*pi/3)+V5s*cos(5*we*t+2*pi/3)+V3s*cos(3*we*t); 
% vcs=V1s*cos(we*t+2*pi/3)+V7s*cos(7*we*t+2*pi/3)+V5s*cos(5*we*t-2*pi/3)+V3s*cos(3*we*t); 

vas=V1s*cos(we*t)        +V29s*cos(29*we*t)        +V25s*cos(25*we*t)        +V23s*cos(23*we*t)        +V19s*cos(19*we*t)        +V17s*cos(17*we*t)        +V13s*cos(13*we*t)        +V11s*cos(11*we*t)        +V7s*cos(7*we*t)        +V5s*cos(5*we*t)        +V3s*cos(3*we*t);
vbs=V1s*cos(we*t-2*pi/3) +V29s*cos(29*we*t+2*pi/3) +V25s*cos(25*we*t-2*pi/3) +V23s*cos(23*we*t+2*pi/3) +V19s*cos(19*we*t-2*pi/3) +V17s*cos(17*we*t+2*pi/3) +V13s*cos(13*we*t-2*pi/3) +V11s*cos(11*we*t+2*pi/3) +V7s*cos(7*we*t-2*pi/3) +V5s*cos(5*we*t+2*pi/3) +V3s*cos(3*we*t-2*pi/3);
vcs=V1s*cos(we*t+2*pi/3) +V29s*cos(29*we*t-2*pi/3) +V25s*cos(25*we*t+2*pi/3) +V23s*cos(23*we*t-2*pi/3) +V19s*cos(19*we*t+2*pi/3) +V17s*cos(17*we*t-2*pi/3) +V13s*cos(13*we*t+2*pi/3) +V11s*cos(11*we*t-2*pi/3) +V7s*cos(7*we*t+2*pi/3) +V5s*cos(5*we*t-2*pi/3) +V3s*cos(3*we*t+2*pi/3);

sf=(we-wr)/we;
wf=sf*we; % equal to wf=we-wr

var=V1r*cos(wf*t)        +V29r*cos(29*wf*t)        +V25r*cos(25*wf*t)        +V23r*cos(23*wf*t)        +V19r*cos(19*wf*t)        +V17r*cos(17*wf*t)        +V13r*cos(13*wf*t)        +V11r*cos(11*wf*t)        +V7r*cos(7*wf*t)        +V5r*cos(5*wf*t); 
vbr=V1r*cos(wf*t-2*pi/3) +V29r*cos(29*wf*t+2*pi/3) +V25r*cos(25*wf*t-2*pi/3) +V23r*cos(23*wf*t+2*pi/3) +V19r*cos(19*wf*t-2*pi/3) +V17r*cos(17*wf*t+2*pi/3) +V13r*cos(13*wf*t-2*pi/3) +V11r*cos(11*wf*t+2*pi/3) +V7r*cos(7*wf*t-2*pi/3) +V5r*cos(5*wf*t+2*pi/3);
vcr=V1r*cos(wf*t+2*pi/3) +V29r*cos(29*wf*t-2*pi/3) +V25r*cos(25*wf*t+2*pi/3) +V23r*cos(23*wf*t-2*pi/3) +V19r*cos(19*wf*t+2*pi/3) +V17r*cos(17*wf*t-2*pi/3) +V13r*cos(13*wf*t+2*pi/3) +V11r*cos(11*wf*t-2*pi/3) +V7r*cos(7*wf*t+2*pi/3) +V5r*cos(5*wf*t-2*pi/3); 

% var=V1r*cos(wf*t)       +V7r*cos(7*wf*t)       +V5r*cos(5*wf*t)       +V3r*cos(3*wf*t); 
% vbr=V1r*cos(wf*t-2*pi/3)+V7r*cos(7*wf*t-2*pi/3)+V5r*cos(5*wf*t+2*pi/3)+V3r*cos(3*wf*t);
% vcr=V1r*cos(wf*t+2*pi/3)+V7r*cos(7*wf*t+2*pi/3)+V5r*cos(5*wf*t-2*pi/3)+V3r*cos(3*wf*t);

vabcs=[vas;vbs;vcs];
vabcr=[var;vbr;vcr];
Vsr=[vabcs;vabcr];

Pabcr=var.*iar+vbr.*ibr+vcr.*icr; % total instantaneous power

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

dwr=P*(Te-Tm)/(2*J); % for motor use (Te-Tm) for generator use (-Te+Tm)
dOr=wr;

% state variables

dX=[dIsr;dwr;dOr];
          
end