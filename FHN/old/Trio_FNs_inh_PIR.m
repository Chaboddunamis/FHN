
function FN2_run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all 
clear all 
%figure('Position',[1 200 1000 1000]); % Specify window size
%global I;   % Need this so that function "f" knows about variable I
eps=0.002;
shift=0.1;
k=10;
k2=3;
I=0.48;
%0.6
g12=.56;
g21=.56;


f = @(t,y) [y(1)-y(1).^3 - y(2) + I-g21*(y(1)+1.5)./(1+exp(-k*(y(5))));eps*(1./(1+exp(-k2*(y(1)+shift)))-y(2));...
            y(3)-y(3).^3 - y(4) + I-g12*(y(3)+1.5)./(1+exp(-k*(y(1))));eps*(1./(1+exp(-k2*(y(3)+shift)))-y(4)); ...
            y(5)-y(5).^3 - y(6) + I-g12*(y(5)+1.5)./(1+exp(-k*(y(3))));eps*(1./(1+exp(-k2*(y(5)+shift)))-y(6)); ];

V1 =[-3:0.01:3];
WN = 1./(1+exp(-k2*(V1+shift)));  % W-nullcline
VN = V1 - V1.^3 + I;   % V-nullcline
VP =V1 - V1.^3 + I -g12*(V1+2);

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  1     %%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,1);
hold off;
plot(WN,V1, 'Color', [.5 .5 .5], 'linewidth', 1);    % plot the W-nullcline
hold on;
plot(VN,V1, 'Color', [.5 .5 .5], 'linewidth', 1);   % plot the V-nullcline
hold on 
axis([-0.1 1.1 -1.3 1.3]);

options = odeset('RelTol',1e-3,'AbsTol',[1e-4],'Events',@events);
[T,Y,tau,Ye,ie] = ode45(f,(0:0.1:20600),[1.5; -0.2; -1.5; 0.5; -1.5; .65],options);

plot(Y(:,2), Y(:,1), 'b', 'linewidth', 3);
hold on
plot(Y(:,4), Y(:,3), 'g', 'linewidth', 3);
hold on 
plot(Y(:,6), Y(:,5), 'r', 'linewidth', 1);
hold on 


t1=tau(ie==1);
t2=tau(ie==2);

Ye1=Ye(:,1);
th1=Ye1(ie==1);

Ye2=Ye(:,3);
th2=Ye2(ie==2);


% t11=t1-circshift(t1,1,2)
% P1=t11(2:end)
 lmin=min([length(t1),length(t2)]);
 
for i=1:lmin-1 
 %   you can find the period of cell 1 like this or see above
P1(i)=t1(i+1)-t1(i);
phaselag1(i)=mod((t2(i)-t1(i))/P1(i),1);
%phaselag2(i)=mod((t3(i)-t1(i))/P1(i),1);
end

title(['Phase plane, I=', num2str(I)], 'fontsize', 16);
xlabel('V(t)'); ylabel('W(t)');
axis([-0.1 1.1 -1.3 1.3]);

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  2     %%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,2); hold off;
plot(T, Y(:,1), 'b', 'linewidth', 1); hold on;
plot (t1,th1,'.', 'MarkerSize',25,'Color',[0./255  81./255  225./255])
hold on

plot(T, Y(:,3), 'g', 'linewidth', 1); hold on;
plot (t2,th2,'.', 'MarkerSize',25,'Color',[0./255  255./255  5./255])
hold on

plot(T, Y(:,5), 'r', 'linewidth', 1); hold on;
hold on

title('Voltage vs time', 'fontsize', 16);
xlabel('time'); ylabel('V(t)');
axis tight;

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  3     %%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,1,3); hold off;

plot(phaselag1,'Color',[2./255  45./255  25./255])
hold on

title('Phase lag vs cycle #', 'fontsize', 16);
xlabel('number'); ylabel('\Delta_{12}');
axis([0 length(t2) -0.05 1.05]);

end

function [value,isterminal,direction] = events(t,y);
th1=y(1);
th2=y(3);
direction= [1,1];
value= [th1,th2];
isterminal=[0,0];
end

