clear all

X0=[0;0;0;0;0;0;0;0];    
T=1/60; N=180; dT=T/N; Ts=0:dT:N*T-dT;      

solution_induction_machine; 
rs=ras; rr=rar;                    
Xls=we*Lls; Xlr=we*Llr;
Lm=3/2*Lms; Xm=we*Lm;

x=size(t);x=x(1); s1=s(x); Oo=Or(x)-wr(x)*t(x); 

% ESTATOR
h=1; sh=s1; 
Vs=V1s;
Vr=V1r*exp(i*Oo);
Vh=[Vs;Vr/sh];
Zh=[rs+j*h*(Xls+Xm)           j*h*Xm
            j*h*Xm rr/sh+j*h*(Xlr+Xm)];    
Ih=inv(Zh)*Vh;
I1s=Ih(1); I1r=Ih(2);
I1st=abs(I1s)*cos(h*we*t+angle(I1s));
I1rt=abs(I1r)*cos(sh*h*we*t+angle(I1r)-Oo);

%  Fun(Hz)      Harm(Hz)       Harm/Fun        Mag 
[s1*we/(2*pi) sh*h*we/(2*pi)  sh*h/s1    abs(I1r)];

h=3; sh=1; 
Vs=V3s;
Vr=0;
I3s=inv(rs+j*h*Xls)*Vs; I3r=inv(rr/sh+j*h*Xlr)*Vr/sh;
I3st=abs(I3s)*cos(h*we*t+angle(I3s));
I3rt=abs(I3r)*cos(sh*h*we*t+angle(I3r));

h=5; sh=s1+wr(x)/we*(1+1/h);
%  Fun(Hz)      Harm(Hz)       Harm/Fun        Mag 
[s1*we/(2*pi) sh*h*we/(2*pi)  sh*h/s1    abs(I1r)];
Vs=V5s;
Vr=0;
Vh=[Vs;Vr/sh];
Zh=[rs+j*h*(Xls+Xm)           j*h*Xm
            j*h*Xm rr/sh+j*h*(Xlr+Xm)];       
Ih=inv(Zh)*Vh;
I5s=Ih(1); I5r=Ih(2);
I5st=abs(I5s)*cos(h*we*t+angle(I5s));
I5rt=abs(I5r)*cos(sh*h*we*t+angle(I5r)+Oo); 

h=7; sh=s1+wr(x)/we*(1-1/h);
%  Fun(Hz)      Harm(Hz)       Harm/Fun        Mag 
[s1*we/(2*pi) sh*h*we/(2*pi)  sh*h/s1    abs(I1r)];
Vs=V7s;
Vr=0;
Vh=[Vs;Vr/sh];
Zh=[rs+j*h*(Xls+Xm)           j*h*Xm
            j*h*Xm rr/sh+j*h*(Xlr+Xm)];       
Ih=inv(Zh)*Vh;
I7s=Ih(1); I7r=Ih(2);
I7st=abs(I7s)*cos(h*we*t+angle(I7s));
I7rt=abs(I7r)*cos(sh*h*we*t+angle(I7r)-Oo); %Use -Oo for positive sequence

h=11; sh=s1+wr(x)/we*(1+1/h);
%  Fun(Hz)      Harm(Hz)       Harm/Fun        Mag 
[s1*we/(2*pi) sh*h*we/(2*pi)  sh*h/s1    abs(I1r)];
Vs=V11s;
Vr=0;
Vh=[Vs;Vr/sh];
Zh=[rs+j*h*(Xls+Xm)           j*h*Xm
            j*h*Xm rr/sh+j*h*(Xlr+Xm)];       
Ih=inv(Zh)*Vh;
I11s=Ih(1); I11r=Ih(2);
I11st=abs(I11s)*cos(h*we*t+angle(I11s));
I11rt=abs(I11r)*cos(sh*h*we*t+angle(I11r)+Oo); 

h=13; sh=s1+wr(x)/we*(1-1/h);
%  Fun(Hz)      Harm(Hz)       Harm/Fun        Mag 
[s1*we/(2*pi) sh*h*we/(2*pi)  sh*h/s1    abs(I1r)];
Vs=V13s;
Vr=0;
Vh=[Vs;Vr/sh];
Zh=[rs+j*h*(Xls+Xm)           j*h*Xm
            j*h*Xm rr/sh+j*h*(Xlr+Xm)];       
Ih=inv(Zh)*Vh;
I13s=Ih(1); I13r=Ih(2);
I13st=abs(I13s)*cos(h*we*t+angle(I13s));
I13rt=abs(I13r)*cos(sh*h*we*t+angle(I13r)-Oo); %Use -Oo for positive sequence

h=17; sh=s1+wr(x)/we*(1+1/h);
%  Fun(Hz)      Harm(Hz)       Harm/Fun        Mag 
[s1*we/(2*pi) sh*h*we/(2*pi)  sh*h/s1    abs(I1r)];
Vs=V17s;
Vr=0;
Vh=[Vs;Vr/sh];
Zh=[rs+j*h*(Xls+Xm)           j*h*Xm
            j*h*Xm rr/sh+j*h*(Xlr+Xm)];       
Ih=inv(Zh)*Vh;
I17s=Ih(1); I17r=Ih(2);
I17st=abs(I17s)*cos(h*we*t+angle(I17s));
I17rt=abs(I17r)*cos(sh*h*we*t+angle(I17r)+Oo); 

h=19; sh=s1+wr(x)/we*(1-1/h);
%  Fun(Hz)      Harm(Hz)       Harm/Fun        Mag 
[s1*we/(2*pi) sh*h*we/(2*pi)  sh*h/s1    abs(I1r)];
Vs=V19s;
Vr=0;
Vh=[Vs;Vr/sh];
Zh=[rs+j*h*(Xls+Xm)           j*h*Xm
            j*h*Xm rr/sh+j*h*(Xlr+Xm)];       
Ih=inv(Zh)*Vh;
I19s=Ih(1); I19r=Ih(2);
I19st=abs(I19s)*cos(h*we*t+angle(I19s));
I19rt=abs(I19r)*cos(sh*h*we*t+angle(I19r)-Oo); %Use -Oo for positive sequence

h=23; sh=s1+wr(x)/we*(1+1/h);
%  Fun(Hz)      Harm(Hz)       Harm/Fun        Mag 
[s1*we/(2*pi) sh*h*we/(2*pi)  sh*h/s1    abs(I1r)];
Vs=V23s;
Vr=0;
Vh=[Vs;Vr/sh];
Zh=[rs+j*h*(Xls+Xm)           j*h*Xm
            j*h*Xm rr/sh+j*h*(Xlr+Xm)];       
Ih=inv(Zh)*Vh;
I23s=Ih(1); I23r=Ih(2);
I23st=abs(I23s)*cos(h*we*t+angle(I23s));
I23rt=abs(I23r)*cos(sh*h*we*t+angle(I23r)+Oo); 

h=25; sh=s1+wr(x)/we*(1-1/h);
%  Fun(Hz)      Harm(Hz)       Harm/Fun        Mag 
[s1*we/(2*pi) sh*h*we/(2*pi)  sh*h/s1    abs(I1r)];
Vs=V25s;
Vr=0;
Vh=[Vs;Vr/sh];
Zh=[rs+j*h*(Xls+Xm)           j*h*Xm
            j*h*Xm rr/sh+j*h*(Xlr+Xm)];       
Ih=inv(Zh)*Vh;
I25s=Ih(1); I25r=Ih(2);
I25st=abs(I25s)*cos(h*we*t+angle(I25s));
I25rt=abs(I25r)*cos(sh*h*we*t+angle(I25r)-Oo); %Use -Oo for positive sequence

h=29; sh=s1+wr(x)/we*(1+1/h);
%  Fun(Hz)      Harm(Hz)       Harm/Fun        Mag 
[s1*we/(2*pi) sh*h*we/(2*pi)  sh*h/s1    abs(I1r)];
Vs=V29s;
Vr=0;
Vh=[Vs;Vr/sh];
Zh=[rs+j*h*(Xls+Xm)           j*h*Xm
            j*h*Xm rr/sh+j*h*(Xlr+Xm)];       
Ih=inv(Zh)*Vh;
I29s=Ih(1); I29r=Ih(2);
I29st=abs(I29s)*cos(h*we*t+angle(I29s));
I29rt=abs(I29r)*cos(sh*h*we*t+angle(I29r)+Oo); 

% Corriente total
% ================================================================
Ist=I1st+I3st+I5st+I7st+I11st+I13st+I17st+I19st+I23st+I25st+I29st;
Irt=I1rt+I3rt+I5rt+I7rt+I11rt+I13rt+I17rt+I19rt+I23rt+I25rt+I29rt;

% ROTOR 
h=1; sh=1; 
%wfr=(1/s1)*we;
wfr=we-wr(x);
Vs=0;%V1s;
Vr=V1r;
Vh=[Vs/sh;Vr];
Zh=[rs/sh+j*h*s1*(Xls+Xm)           j*h*s1*Xm
            j*h*s1*Xm         rr+j*h*s1*(Xlr+Xm)];     
Ih=inv(Zh)*Vh;
I1sx=Ih(1); I1rx=Ih(2);
I1stx=abs(I1sx)*cos(h*we*t+angle(I1sx)+Oo+.1);
I1rtx=abs(I1rx)*cos(s1*h*we*t+angle(I1rx)*.11);

h=5; sh=(-h*wfr+wr(x))/(-wfr);
Vs=0;
Vr=V5r;
Vh=[Vs/sh;Vr];
Zh=[rs/sh+j*h*s1*(Xls+Xm)           j*h*s1*Xm
            j*h*s1*Xm         rr+j*h*s1*(Xlr+Xm)];     
Ih=inv(Zh)*Vh;
I5sx=Ih(1); I5rx=Ih(2);
I5stx=abs(I5sx)*cos(sh*h*we*t+angle(I5sx)-Oo);
I5rtx=abs(I5rx)*cos(s1*h*we*t+angle(I5rx));

h=7; sh=s1+wr(x)/we*(1-h);
Vs=0;
Vr=V7r;
Vh=[Vs/sh;Vr];
Zh=[rs/sh+j*h*s1*(Xls+Xm)           j*h*s1*Xm
            j*h*s1*Xm         rr+j*h*s1*(Xlr+Xm)];      
Ih=inv(Zh)*Vh;
I7sx=Ih(1); I7rx=Ih(2);
I7stx=abs(I7sx)*cos(sh*h*we*t+angle(I7sx)+Oo);
I7rtx=abs(I7rx)*cos(s1*h*we*t+angle(I7rx));

h=11; sh=(-h*wfr+wr(x))/(-wfr);
Vs=0;
Vr=V11r;
Vh=[Vs/sh;Vr];
Zh=[rs/sh+j*h*s1*(Xls+Xm)           j*h*s1*Xm
            j*h*s1*Xm         rr+j*h*s1*(Xlr+Xm)];     
Ih=inv(Zh)*Vh;
I11sx=Ih(1); I11rx=Ih(2);
I11stx=abs(I11sx)*cos(sh*h*we*t+angle(I11sx)-Oo);
I11rtx=abs(I11rx)*cos(s1*h*we*t+angle(I11rx));

h=13; sh=s1+wr(x)/we*(1-h);
Vs=0;
Vr=V13r;
Vh=[Vs/sh;Vr];
Zh=[rs/sh+j*h*s1*(Xls+Xm)           j*h*s1*Xm
            j*h*s1*Xm         rr+j*h*s1*(Xlr+Xm)];      
Ih=inv(Zh)*Vh;
I13sx=Ih(1); I13rx=Ih(2);
I13stx=abs(I13sx)*cos(sh*h*we*t+angle(I13sx)+Oo);
I13rtx=abs(I13rx)*cos(s1*h*we*t+angle(I13rx));

h=17; sh=(-h*wfr+wr(x))/(-wfr);
Vs=0;
Vr=V17r;
Vh=[Vs/sh;Vr];
Zh=[rs/sh+j*h*s1*(Xls+Xm)           j*h*s1*Xm
            j*h*s1*Xm         rr+j*h*s1*(Xlr+Xm)];     
Ih=inv(Zh)*Vh;
I17sx=Ih(1); I17rx=Ih(2);
I17stx=abs(I17sx)*cos(sh*h*we*t+angle(I17sx)-Oo);
I17rtx=abs(I17rx)*cos(s1*h*we*t+angle(I17rx));

h=19; sh=s1+wr(x)/we*(1-h);
Vs=0;
Vr=V19r;
Vh=[Vs/sh;Vr];
Zh=[rs/sh+j*h*s1*(Xls+Xm)           j*h*s1*Xm
            j*h*s1*Xm         rr+j*h*s1*(Xlr+Xm)];      
Ih=inv(Zh)*Vh;
I19sx=Ih(1); I19rx=Ih(2);
I19stx=abs(I19sx)*cos(sh*h*we*t+angle(I19sx)+Oo);
I19rtx=abs(I19rx)*cos(s1*h*we*t+angle(I19rx));

h=23; sh=(-h*wfr+wr(x))/(-wfr);
Vs=0;
Vr=V23r;
Vh=[Vs/sh;Vr];
Zh=[rs/sh+j*h*s1*(Xls+Xm)           j*h*s1*Xm
            j*h*s1*Xm         rr+j*h*s1*(Xlr+Xm)];     
Ih=inv(Zh)*Vh;
I23sx=Ih(1); I23rx=Ih(2);
I23stx=abs(I23sx)*cos(sh*h*we*t+angle(I23sx)-Oo);
I23rtx=abs(I23rx)*cos(s1*h*we*t+angle(I23rx));

h=25; sh=s1+wr(x)/we*(1-h);
Vs=0;
Vr=V25r;
Vh=[Vs/sh;Vr];
Zh=[rs/sh+j*h*s1*(Xls+Xm)           j*h*s1*Xm
            j*h*s1*Xm         rr+j*h*s1*(Xlr+Xm)];      
Ih=inv(Zh)*Vh;
I25sx=Ih(1); I25rx=Ih(2);
I25stx=abs(I25sx)*cos(sh*h*we*t+angle(I25sx)+Oo);
I25rtx=abs(I25rx)*cos(s1*h*we*t+angle(I25rx));

h=29; sh=(-h*wfr+wr(x))/(-wfr);
Vs=0;
Vr=V29r;
Vh=[Vs/sh;Vr];
Zh=[rs/sh+j*h*s1*(Xls+Xm)           j*h*s1*Xm
            j*h*s1*Xm         rr+j*h*s1*(Xlr+Xm)];     
Ih=inv(Zh)*Vh;
I29sx=Ih(1); I29rx=Ih(2);
I29stx=abs(I29sx)*cos(sh*h*we*t+angle(I29sx)-Oo);
I29rtx=abs(I29rx)*cos(s1*h*we*t+angle(I29rx));

% ESTATOR + ROTOR
Ist=Ist+I1stx+I5stx+I7stx+I11stx+I13stx+I17stx+I19stx+I23stx+I25stx+I29stx;
Irt=Irt+I1rtx+I5rtx+I7rtx+I11rtx+I13rtx+I17rtx+I19rtx+I23rtx+I25rtx+I29rtx;

% GRAFICACION
figure(4);
subplot(211);
plot(t,Ist,t,ias);%axis([0 3 -12 20]);
%legend('Steady-State','Dynamic');
xlabel('Time [sec]');ylabel('Stator Current');
subplot(212);
plot(t,Irt,t,iar);%axis([0 3 -8 10]);
%legend('Steady-State','Dynamic','Measurement');
xlabel('Time [sec]');ylabel('Rotor Current');
p=size(t);p=x(1); Istc=Ist(p-N+1:p); iasc=ias(p-N+1:p); tx=0:dT:N*dT-dT;
F=tdf(Istc,N);[Ish,h,H]=ordena(F,N); F=tdf(iasc,N);[Iash,h,H]=ordena(F,N);

p=size(t);p=x(1); Irtc=Irt(p-N+1:p); iarc=iar(p-N+1:p); tx=0:dT:N*dT-dT;
F=tdf(Irtc,N);[Irh,h,H]=ordena(F,N); F=tdf(iarc,N);[Iarh,h,H]=ordena(F,N);

figure(5)
subplot(211);
plot(tx,Istc,'-.k',tx,iasc,'k','linewidth',2)
legend('Steady-State','Dynamic');
xlabel('Time [sec]');ylabel('Stator Current');
subplot(212);
bar(h*fe,abs([transpose(Ish) transpose(Iash)]),1);%axis([0 500 0 2.3]);
legend('Steady-State','Dynamic');
xlabel('Frequency [Hz]');ylabel('I_s [A]');
frx=fr(p);Nr=round(fe*N/frx);Irtc=Irt(p-Nr+1:p);iarc=iar(p-Nr+1:p);tx=0:dT:Nr*dT-dT; 
F=tdf(Irtc,Nr);[Irh,h,H]=ordena(F,Nr); F=tdf(iarc,Nr);[Iarh,h,H]=ordena(F,Nr);

figure(6)
subplot(211);
plot(tx,Irtc,'-.k',tx,iarc,'k');%axis([0 0.45 -2 3]);
legend('Steady-State','Dynamic');
xlabel('Time [sec]');ylabel('Rotor Current');
subplot(212);
bar(h*frx,abs([transpose(Irh) transpose(Iarh)]),1)
legend('Steady-State','Dynamic');
xlabel('Frequency [Hz]');ylabel('I_r [A]');