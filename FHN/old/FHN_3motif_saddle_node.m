% 
% dv/dt = -v*(v-a)*(v-1)-w+I0
% dw/dt = eps*(v-d*w)
% 
% This code also plots nullclines of the FHN eqns.

function FHN
clear all
%clf 
% dimensionless model of neuronal excitability
% v=Y(1); w=Y(2);
v0=[-1; -1; 1];w0=[-.08 ;-.232; .123];
a=0.1;d=1;eps=0.05;I0=0.0045;
Y0=[v0;w0];
t=0:1:3990;
options=odeset('RelTol',1.e-5);
[T, Y]=ode45(@dydt_FHN,t,Y0,options,a,eps,d,I0);

len=length(T);
l1=1;
figure(1);clf;

subplot(3,1,1)
plot(T(l1:len),Y(l1:len,1),T(l1:len),Y(l1:len,2),T(l1:len),Y(l1:len,3)); % time courses of V and w
hold on;


[psor.a,lsor.a]=findpeaks(Y(:,1),'MinPeakHeight',0.5);
for i=1:size(lsor.a)
    j=lsor.a(i);
    lcs.a(i)=T(j);
end
%plot(lcs.a,psor.a,'*b');


[psor.b,lsor.b]=findpeaks(Y(:,2),'MinPeakHeight',0.5);
for i=1:size(lsor.b)
    j=lsor.b(i);
    lcs.b(i)=T(j);
end
%plot(lcs.b,psor.b,'*g');

[psor.c,lsor.c]=findpeaks(Y(:,3),'MinPeakHeight',0.5);
for i=1:size(lsor.c)
    j=lsor.c(i);
    lcs.c(i)=T(j);
end
%plot(lcs.c,psor.c,'*r');
length(lcs.a)
length(lcs.b)
length(lcs.c) 

i=min([length(lcs.a),length(lcs.b),length(lcs.c)]);
for j=1:i-1
    tn.a(j)=lcs.a(j+1)-lcs.a(j);
    tn.b(j)=lcs.b(j+1)-lcs.b(j);
    tn.c(j)=lcs.c(j+1)-lcs.c(j);
end
for j=1:i-1
tau.ab(j)=mod((lcs.b(j)-lcs.a(j))/tn.a(i-1),1);
tau.bc(j)=mod((lcs.c(j)-lcs.a(j))/tn.a(i-1),1);
end
% 
legend('v1(t)','v2(t)','v3(t)');
xlabel('Time'); ylabel('v1, v2, v3');
vpts=(-1.5:.05:1.5);
%figure(2);clf;
subplot(3,1,2)
hold on;

plot(Y(:,4),Y(:,1)); % V-w phase plane 
plot(Y(:,5),Y(:,2));
plot(Y(:,6),Y(:,3),'Linewidth',1);
%determine and plot the v,w-nullclines
%options=optimset; % sets options in fzero to default values
%for k=1:61
%vnullpts(k)=fzero(@vrhs_FHN,[-10 10],options,vpts(k),a,I0);  
%end
vnullpts=-vpts.*(vpts-a).*(vpts-1)+I0;
vpts1 =[-3:.1:3];
wnullpts=1./(1+exp(-10*(vpts1-0.75)));
plot(vnullpts,vpts,'black',wnullpts,vpts1,'green');
xlabel('w'); ylabel('v');
axis([-.05 .25 -1 1.5]);
axis on;
subplot(3,1,3)
hold on;

plot(tau.ab,'color','blue');
hold on
plot(tau.bc,'color','magenta');
hold on
title('Phase lags','FontSize',18)
 % legend('\Delta\phi_{1,2}^n','\Delta \phi_{2,3}^n','Location','northeast',24)
    xlabel('cycle number n','FontSize',18)
  ylabel('\Delta \phi^n','FontSize',18)
  XY=ones(i-1)/3;
plot(XY,'-.b');
axis([0 i-1 0 1]);
%plot(tau.bc,'color','red');

figure(2)
plot3(Y(l1:len,4),Y(l1:len,5),Y(l1:len,6)); 

end


function val=vrhs_FHN(w,v,a,I0)
	val=-v*(v-a)*(v-1)-w+I0;
end


function dY=dydt_FHN(t,Y,a,eps,d,I0)
v(1)=Y(1);
v(2)=Y(2);
v(3)=Y(3);
w(1)=Y(4);
w(2)=Y(5);
w(3)=Y(6);
lambda=-11;
thetas=0.5;
VS=-0.7;
g(1)=.005;
g(2)=.0051;
g(3)=.005;
dY=zeros(6,1);
dv=zeros(3,1);
dw=zeros(3,1);
S(1)=0.03*sin(1*t)  +1/(1+exp(lambda*(v(2)-thetas)))+1/(1+exp(lambda*(v(3)-thetas)));
S(2)=0.05*sin(1.2*t)+1/(1+exp(lambda*(v(1)-thetas)))+1/(1+exp(lambda*(v(3)-thetas)));
S(3)=0.075*sin(4.2*t) +1/(1+exp(lambda*(v(1)-thetas)))+1/(1+exp(lambda*(v(2)-thetas)));
dv=-v.*(v-a).*(v-1)-w+I0-g.*(v-VS).*S;
dw=eps*(-w+1./(1+exp(-10*(v-0.75))));
dY(1)=dv(1);
dY(2)=dv(2);
dY(3)=dv(3);
dY(4)=dw(1);
dY(5)=dw(2);
dY(6)=dw(3);

end

