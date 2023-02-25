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
v0=[-1.5; -1.5; -1.5];w0=[0.2 ;0.21; 0.2117];
a=0.1;d=1;eps=0.003;I0=.06;
Y0=[v0;w0];
t=0:0.1:25000;
options=odeset('RelTol',1.e-5);
[T, Y]=ode45(@dydt_FHN,t,Y0,options,a,eps,d,I0);

len=length(T);
l1=len-10000;
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
% plot(lcs.c,psor.c,'*r');
length(lcs.a)
length(lcs.b)
length(lcs.c) 

i=min([length(lcs.a),length(lcs.b),length(lcs.c)])
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
plot(Y(:,6),Y(:,3));
%determine and plot the v,w-nullclines
%options=optimset; % sets options in fzero to default values
%for k=1:61
%vnullpts(k)=fzero(@vrhs_FHN,[-10 10],options,vpts(k),a,I0);  
%end
vnullpts=-vpts.*(vpts-a).*(vpts-1)+I0;
wnullpts=vpts/d;
plot(vnullpts,vpts,'black',wnullpts,vpts,'green');
xlabel('v'); ylabel('w');
axis([0.0 0.25 -.7 1.2]);
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
lambda=-15;
thetas=0.3;
VS=-0.7;
g(1)=.0005;
g(2)=.0005;
g(3)=.0005;
dY=zeros(6,1);
dv=zeros(3,1);
dw=zeros(3,1);
S(1)=1/(1+exp(lambda*(v(2)-thetas)))+1/(1+exp(lambda*(v(3)-thetas)));
S(2)=1/(1+exp(lambda*(v(1)-thetas)))+1/(1+exp(lambda*(v(3)-thetas)));
S(3)=1/(1+exp(lambda*(v(1)-thetas)))+1/(1+exp(lambda*(v(2)-thetas)));
dv=-v.*(v-a).*(v-1)-w+I0*1/(1+exp(20-t)/.2)-g.*(v-VS).*S;
dw=eps*(v-d*w);
dY(1)=dv(1);
dY(2)=dv(2);
dY(3)=dv(3);
dY(4)=dw(1);
dY(5)=dw(2);
dY(6)=dw(3);

%dY=[dv;    dw];


% dY(1)=-v1*(v1-a)*(v1-1)-w1+I0*1/(1+exp(20-t)/.2);
% dY(2)=-v2*(v2-a)*(v2-1)-w2+I0*1/(1+exp(20-t)/.2);
% dY(3)=-v3*(v3-a)*(v3-1)-w3+I0*1/(1+exp(20-t)/.2);
% dY(4)=eps*(v1-d*w1);
% dY(5)=eps*(v1-d*w1);
% dY(6)=eps*(v3-d*w3);
end

