
clear all 
tic
tmax=2000; % max time in seconds
step=0.2;       % time step in msec

%-----Cellular parameters ---------------
eps=0.02;
Vshift=-0.3;
slope=.8;

% Constant stimuli if any
%gsyn=0;
Iext=0.04;
t1=950; % this time is in msec
t2=1200; % this time is in msec
 
time=zeros(tmax/step,1);vv1=time;vv2=time;Caa1=time;Caa2=time;

V1= -2; V2= -2; Ca1=.9; Ca2=0.1; 
tt=0;
i=0;
while (tt < tmax)
%pulse
%Iapp=-gsyn*(V1+2)*heaviside(tt-t1)*heaviside(t2-tt)
%Iapp=Iext*heaviside(tt-t1)*heaviside(t2-tt);
 if tt>t1 && tt<t2
 Iapp=Iext;
 else
 Iapp=0;
 end
    
V1 =V1 +step*( V1-V1.^3 - Ca1 + Iapp );
V2 =V2 +step*( V2-V2.^3 - Ca2 + Iapp );
Ca1=Ca1+step*(eps*(V1-Vshift-slope*Ca1 ));
Ca2=Ca2+step*(eps*(V2-Vshift-slope*Ca2));

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

%ss1(i)=s1;
%ss2(i)=s2;

end  
 
toc     
figure(1)
clf
subplot(3,1,1)
plot(time,vv1,'Color',[0 0  .7],'LineWidth',1.5)
hold on
xlim([0 tmax]) 
ylim([-2 2])
xlabel('Time'),ylabel('Voltage')
 
subplot(3,1,2)
plot(time,vv2,'Color',[0 .7 0],'LineWidth',1.5)
hold on
xlabel('Time'),ylabel('Voltage')
xlim([0 tmax]) 
ylim([-2 2])

subplot(3,1,3)
plot(time,Caa1,'Color',[0 0  .7],'LineWidth',1.5)
hold on
plot(time,Caa2,'Color',[0 .7  0],'LineWidth',1.5)
hold on
%ylim([0.5 1.05])
xlim([0 tmax])
xlabel('Time'),ylabel('[Ca]','Fontsize', 16)
 

 
%   subplot(4,1,3)
%  % plot(time,10*mm1.*ss1,'blue')
%   hold on
%   plot(time,ss1,'Color',[0.2 0.2 0.99])
%   hold on
%   plot(time,ss2,'green')
%   hold on
%   %ylim([0 .02])
%   xlim([0 time(end)])
%   

figure(2)
clf
plot(Caa1,vv1,'Color',[0 0  1],'LineWidth',1.5)
hold on
  
V = [-2:0.1:2];
VN = V - V.^3 + Iext;   % V-nullcline
plot (VN,V,'--','Color',[0.5 0.5 0.5])
hold on
VN = V - V.^3 + 0*Iext;   % V-nullcline
plot (VN,V,'Color',[0.5 0.5 0.5])
hold on
CaN=(V-Vshift)/slope; % Ca-nullcline
plot (CaN,V,'Color',[0.5 0.2 0.2])
hold on 
xlim([-1 1]) 
ylim([-2 2])
 xlabel('[Ca]-variable','Fontsize', 16),ylabel('Voltage','Fontsize', 16)