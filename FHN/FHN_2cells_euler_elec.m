
clear all 
tic
tmax=2000; % max time in seconds
step=0.2;       % time step in msec

%-----Cellular parameters ---------------
eps=0.02;
Vshift=-0.2;
alpha=.1;
beta=0.01;


% Constant stimuli if any
gsyn=0.005;

Iext=0.2;
Iapp=Iext;
t1=100; % this time is in msec
t2=2000; % this time is in msec
 

time=zeros(tmax/step,1);vv1=time;vv2=time;Caa1=time;Caa2=time;

V1= -2; V2= 2; Ca1=-.8; Ca2=0.81; 

tt=0;
i=0;
while (tt < tmax)
    
V1 =V1 +step*( V1-V1.^3 - Ca1 + Iapp - gsyn*(V1-V2));
V2 =V2 +step*( V2-V2.^3 - Ca2 + Iapp  -gsyn*(V2-V1));
Ca1=Ca1+step*(eps*(V1-Vshift-Ca1 ));
Ca2=Ca2+step*(eps*(V2-Vshift-Ca2));
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
 

 
% subplot(4,1,4)
%    plot(time,ss1,'Color',[0.0 0.0 0.7],'LineWidth',1.5)
%    hold on
%    plot(time,ss2,'Color',[0.0 0.6 0.0],'LineWidth',1.5)
%    hold on
% %   %ylim([0 .02])
%    xlim([0 tmax])
%   

figure(2)
clf
plot(Caa1,vv1,'Color',[0 0  1],'LineWidth',1.5)
hold on
  
V = [-2:0.01:2];
VN =V - V.^3  +Iext;   % V-nullcline
plot (VN,V,'.','Color',[0.5 0.5 0.5])
hold on
CaN=(V-Vshift); % Ca-nullcline
plot (CaN,V,'Color',[0.8 0.2 0.2])
hold on 
xlim([-1 1]) 
ylim([-2 2])
 xlabel('[Ca]-variable','Fontsize', 16),ylabel('Voltage','Fontsize', 16)