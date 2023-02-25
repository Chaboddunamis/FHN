%
%                     Fitzhugh-Nagumo Model
%          Subcritical and supercritical Hopf Bifurcations
%                
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf 
close all 
figure('Position',[1 200 1100 1000]); % Specify window size
%global I;   % Need this so that function "f" knows about variable I
eps=0.1;
shift=.1;
k=5;
I=0.1
%for I = [1.4:-0.01:0.6];  % Loop over different values of applied current I

    for I = [.68:0.0001:.7];
f = @(t,y) [ y(1) - y(1).^3/3 - y(2) + I; eps*(1/(1+exp(-k*(y(1)+shift)))- 0.5*y(2)) ];
g = @(y) f(0,y);

                      % Eigenvalues of Jacobian

V = [-3:0.02:3];
WN = 1./(1+exp(-k*(V+shift)))/0.5;  % W-nullcline
VN = V - V.^3/3 + I;   % V-nullcline

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  1     %%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1);
hold off;
plot(WN,V, 'r-.', 'linewidth', 1);    % plot the W-nullcline
hold on;
plot(VN,V, 'k', 'linewidth', 1);   % plot the V-nullcline
hold on 


[T Y] = ode45(f, [0 2200], [-1, 2]);
plot(Y(:,2), Y(:,1), 'b-', 'linewidth', 1.5);

opts = optimoptions('fminunc','Display','none','Algorithm','quasi-newton');
fp = fsolve(g, [Y(end,1), Y(end,2)],opts); 
%fp = fsolve(g, [0, -2],opts); % Find the fixed point

Vss = fp(1); Wss = fp(2);   % Get the steady-state V and W values from "fp"
J = [ [1 - Vss^2, -1]; [eps*(exp(-k*(Vss+shift)))*(-k)./(1.+exp(-k*(Vss+shift)))^2,  -0.5*eps ]];  % The Jacobian
Lambda = eig(J);    
plot(Wss,Vss, '*', 'color','g');
hold on 
plot(Y(end,2),Y(end,1), '*', 'color','b');
hold on 

title(['Phase plane, I=', num2str(I)], 'fontsize', 16);
xlabel('V(t)'); ylabel('W(t)');
axis([-0.5 2.5 -2.5 2.5]);

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  2     %%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2); hold off;
plot(T, Y(:,1), 'b-', 'linewidth', 1.5); hold on;
plot([0 2000],[Vss Vss],'b:', 'linewidth', 1.5);
title('Voltage vs time', 'fontsize', 16);
xlabel('time'); ylabel('V(t)');
axis tight;

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  3     %%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3); hold on;

%if real(Lambda(1)) > 0
    color1 = [1 0 1];   % Red color  ([R G B] values)
%else 
    color2 = [0 1 1];   % Blue color ([R G B] values)
%end;
plot(real(Lambda(1)),imag(Lambda(1)), '*', 'color', color1);
plot(real(Lambda(2)),imag(Lambda(2)), '*', 'color', color2);
title(['Real(\lambda)=',num2str(real(Lambda(2)))], 'fontsize', 12);
axis tight;

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  4     %%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4); hold on;
[T Y] = ode45(f, [0 200], Y(length(T),:) );  % Run some more in order;  % to settle at equilibrium
plot(I, min(Y(:,1)), 'g.');       % minimum of V(t): 
plot(I, max(Y(:,1)), 'b.');       % maximum of V(t): 
plot(I,Vss,'color',color2);    % Fixed point: a red or blue point
title('Bifurcation diagram (V_{min}, V_{max}, V*)', 'fontsize', 12);
xlabel('I'); ylabel('V equilibrium');
axis tight;
drawnow;

end;



