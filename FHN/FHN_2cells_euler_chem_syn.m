
clear all 
tic
tmax=2000; % max time in seconds
step=0.2;       % time step in msec

%-----Cellular parameters ---------------
eps=0.02;
Vshift=-0.2;
alpha=1;
beta=0.01;


% Constant stimuli if any
gsyn12=0.02;
gsyn21=0.02;

Erev=-2;

Iext=-0.1;
t1=600; % this time is in msec
t2=650; % this time is in msec
 

time=zeros(tmax/step,1);vv1=time;vv2=time;Caa1=time;Caa2=time;
ss1=time;ss2=time;

V1= 2; V2= -2; Ca1=.61; Ca2=-0.6; s1=0; s2=0;

tt=0;
i=0;
while (tt < tmax)
%pulse
 if tt>t1 && tt<t2
 Iapp=Iext;
 else
      Iapp=0.;
 end
    
V1 =V1 +step*( V1-V1.^3 - Ca1 + Iapp -gsyn21*s2*(V1-Erev));
V2 =V2 +step*( V2-V2.^3 - Ca2        - gsyn12*s1*(V2-Erev));
Ca1=Ca1+step*(eps*(V1-Vshift-Ca1 ));
Ca2=Ca2+step*(eps*(V2-Vshift-Ca2));
s1=s1  +step*(alpha*(1-s1)/(1+exp(-50*(V1))))-beta*s1;
s2=s2  +step*(alpha*(1-s2)/(1+exp(-50*(V2))))-beta*s2;
tt=tt+step;
i=i+1;
time(i)=tt;
vv1(i)=V1;
vv2(i)=V2;
Caa1(i)=Ca1;
Caa2(i)=Ca2;

%plot (tt,V1,'.','Color','blue');
%hold on
%drawnow; 

ss1(i)=s1;
ss2(i)=s2;

end  
 
toc     
figure(1)
clf
subplot(4,1,1)
plot(time,vv1,'Color',[0 0  .7],'LineWidth',1.5)
hold on
xlim([0 tmax]) 
ylim([-2 2])
xlabel('Time'),ylabel('Voltage')
 
subplot(4,1,2)
plot(time,vv2,'Color',[0 .7 0],'LineWidth',1.5)
hold on
xlabel('Time'),ylabel('Voltage')
xlim([0 tmax]) 
ylim([-2 2])

subplot(4,1,3)
plot(time,Caa1,'Color',[0 0  .7],'LineWidth',1.5)
hold on
plot(time,Caa2,'Color',[0 .7  0],'LineWidth',1.5)
hold on
%ylim([0.5 1.05])
xlim([0 tmax])
xlabel('Time'),ylabel('[Ca]','Fontsize', 16)
 

 
subplot(4,1,4)
   plot(time,ss1,'Color',[0.0 0.0 0.7],'LineWidth',1.5)
   hold on
   plot(time,ss2,'Color',[0.0 0.6 0.0],'LineWidth',1.5)
   hold on
   ylim([0 1.1])
   xlim([0 tmax])
%   

figure(2)
clf
plot(Caa1,vv1,'Color',[0 0  1],'LineWidth',1.5)
hold on
  
V = [-2:0.1:2];
VN =V - V.^3 + -gsyn12*(V1-Erev);   % V-nullcline
plot (VN,V,'Color',[0.5 0.5 0.5])
hold on
VN = V - V.^3 + 0*Iext;   % V-nullcline
plot (VN,V,'Color',[0.2 0.2 0.2])
hold on
CaN=(V-Vshift); % Ca-nullcline
plot (CaN,V,'Color',[0.8 0.2 0.2])
hold on 
xlim([-1 1]) 
ylim([-2 2])
 xlabel('[Ca]-variable','Fontsize', 16),ylabel('Voltage','Fontsize', 16)