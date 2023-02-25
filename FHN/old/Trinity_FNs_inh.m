
function FN2_run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
clear all 
figure('Position',[1 200 1000 1000]); % Specify window size
%global I;   % Need this so that function "f" knows about variable I
eps=0.005;
shift=0;
k=10;
I=0.4;
%0.6
g=0.03;

f = @(t,y) [y(1)-y(1).^3 - y(2) + I-g*(y(1)+2).*(1./(1+exp(-k*(y(3))))+1./(1+exp(-k*(y(5))))); eps*(1./(1+exp(-k*(y(1)+shift)))-y(2)); ....
            y(3)-y(3).^3 - y(4) + I-g*(y(3)+2).*(1./(1+exp(-k*(y(1))))+1./(1+exp(-k*(y(5))))); eps*(1./(1+exp(-k*(y(3)+shift)))-y(4)); ...
            y(5)-y(5).^3 - y(6) + I-g*(y(3)+2).*(1./(1+exp(-k*(y(1))))+1./(1+exp(-k*(y(3))))); eps*(1./(1+exp(-k*(y(5)+shift)))-y(6)) ];

V1 =[-3:0.01:3];
WN = 1./(1+exp(-k*(V1+shift)));  % W-nullcline
VN = V1 - V1.^3 + I;   % V-nullcline
VP =V1 - V1.^3 + I -g*(V1+2);

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  1     %%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,1);
hold off;
plot(WN,V1, 'Color', [.5 .5 .5], 'linewidth', 1);    % plot the W-nullcline
hold on;
plot(VN,V1, 'Color', [.5 .5 .5], 'linewidth', 1);   % plot the V-nullcline
hold on 
plot(VP,V1, 'Color', [.5 .0 .0], 'linewidth', 1);   % plot the V-nullcline
hold on 


options = odeset('RelTol',1e-4,'AbsTol',[1e-4],'Events',@events);
[T,Y,tau,Ye,ie] = ode45(f,(0:0.1:56000),[0; 1.0222; 0; 1.022; 0; 1.0221],options);


plot(Y(:,2), Y(:,1), 'b', 'linewidth', 1);
hold on
plot(Y(:,4), Y(:,3), 'g', 'linewidth', 1);
hold on 
plot(Y(:,6), Y(:,5), 'red', 'linewidth', 1);
hold on 

t1=tau(ie==1);
t2=tau(ie==2);
t3=tau(ie==3);

Ye1=Ye(:,1);
th1=Ye1(ie==1);

Ye2=Ye(:,3);
th2=Ye2(ie==2);

Ye3=Ye(:,5);
th3=Ye3(ie==3);



% t11=t1-circshift(t1,1,2)
% P1=t11(2:end)
 lmin=min([length(t1),length(t2),length(t3)]);
 
for i=1:lmin-1 
 %   you can find the period of cell 1 like this or see above
P1(i)=t1(i+1)-t1(i);
phaselag1(i)=mod((t2(i)-t1(i))/P1(i),1); %#ok<AGROW>
phaselag2(i)=mod((t3(i)-t1(i))/P1(i),1);
end

title(['Phase plane, I=', num2str(I)], 'fontsize', 16);
xlabel('V(t)'); ylabel('W(t)');
axis([-0.1 1.1 -1.3 1.3]);

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  2     %%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,2); hold off;
plot(T, Y(:,1), 'b', 'linewidth', 1); hold on;
plot (t1,0,'.', 'MarkerSize',25,'Color',[0./255  81./255  225./255])
hold on

plot(T, Y(:,3), 'green', 'linewidth', 1); hold on;
plot (t2,0,'.', 'MarkerSize',25,'Color',[0./255  255./255  5./255])
hold on

plot(T, Y(:,5), 'red', 'linewidth', 1); hold on;
plot (t3,0,'.', 'MarkerSize',25,'Color',[0./255  255./255  5./255])
hold on

title('Voltage vs time', 'fontsize', 16);
xlabel('time'); ylabel('V(t)');
axis tight;

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  3     %%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,1,3); hold off;

plot(phaselag1,'Color',[2./255  245./255  25./255])
hold on
plot(phaselag2,'Color',[225./255  2./255  25./255])
hold on

title('Phase lags vs cycle #', 'fontsize', 16);
xlabel('number'); ylabel('\Delta_{12} and \Delta_{13}');
axis([0 length(t2) -0.05 1.05]);

end

function [value,isterminal,direction] = events(t,y);
th1=y(1);
th2=y(3);
th3=y(5);
direction= [1,1,1];
value= [th1,th2,th3];
isterminal=[0,0,0];
end

