% solution induction_machine

% zero initial conditions, simulation time 0.6 sec.

[t,X]=ode23('induction_machine',Ts,X0);
induction_machine_23HP;

% note: the voltage magnitudes are in the data file of the machine

vas=V1s*cos(we*t)       +V7s*cos(7*we*t)        +V5s*cos(5*we*t)        +V3s*cos(3*we*t); 
vbs=V1s*cos(we*t-2*pi/3)+V7s*cos(7*we*t-2*pi/3) +V5s*cos(5*we*t+2*pi/3) +V3s*cos(3*we*t); 
vcs=V1s*cos(we*t+2*pi/3)+V7s*cos(7*we*t+2*pi/3) +V5s*cos(5*we*t-2*pi/3) +V3s*cos(3*we*t);

ias=X(:,1); ibs=X(:,2); ics=X(:,3);
iar=X(:,4); ibr=X(:,5); icr=X(:,6);
wr =X(:,7); Or =X(:,8);

i1=ias.*(iar-ibr/2-icr/2)+ibs.*(ibr-iar/2-icr/2)+ics.*(icr-ibr/2-iar/2);
i2=ias.*(ibr-icr)+ibs.*(icr-iar)+ics.*(iar-ibr);
Te=-(P/2)*Lms*(i1.*sin(Or)+sqrt(3)/2*i2.*cos(Or));

wrm=60/(pi*P)*wr;  % motor speed in r.p.m.
s=(we-wr)/we;      % slip
fr=s*fe;           % frecuency of the currents in the rotor
Pabcs=vas.*ias+vbs.*ibs+vcs.*ics; % total instantaneous power

figure(1);
subplot(311);plot(t,ias);ylabel('ias [Amp]');xlabel('time [sec]');
subplot(312);plot(t,ibs);ylabel('ibs [Amp]');xlabel('time [sec]');
subplot(313);plot(t,ics);ylabel('ics [Amp]');xlabel('time [sec]');

figure(2);
subplot(311);plot(t,iar);ylabel('iar [Amp]');xlabel('time [sec]');
subplot(312);plot(t,ibr);ylabel('ibr [Amp]');xlabel('time [sec]');
subplot(313);plot(t,icr);ylabel('icr [Amp]');xlabel('time [sec]');

figure(3);
subplot(311);plot(t,Te);ylabel('Electric torque Te [N.m]');xlabel('time [sec]');
subplot(312);plot(t,wrm);ylabel('motor speed wrm [r.p.m]');xlabel('time [sec]');
subplot(313);plot(wrm,Te);ylabel('Torque Te');xlabel('speed wrm');

