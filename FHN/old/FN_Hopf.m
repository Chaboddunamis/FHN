%
%                     Fitzhugh-Nagumo Model
%          Subcritical and supercritical Hopf Bifurcations
%                
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf 
close all 
figure('Position',[1 200 1200 1000]); % Specify window size
%global I;   % Need this so that function "f" knows about variable I
eps=0.1;
for I = [0.38:0.01:.7];  % Loop over different values of applied current I

f = @(t,y) [ (y(1) - y(1).^3/3 - y(2) + I); eps*(y(1) + 0.7 - y(2)) ];
f1 = @(t,y) [ -(y(1) - y(1).^3/3 - y(2) + I); -eps*(y(1) + 0.7 - y(2)) ];
g = @(y) f(0,y);

fp = fsolve(g,[0 0]);       % Find the fixed point
Vss = fp(1); Wss = fp(2);   % Get the steady-state V and W values from "fp"
J = [ [1 - Vss^2, -1]; [eps, -eps]];  % The Jacobian
Lambda = eig(J);                          % Eigenvalues of Jacobian

V = [-3:0.05:3];
WN = (V + 0.7) / 1;  % W-nullcline
VN = V - V.^3/3 + I;   % V-nullcline

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  1     %%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1);
hold off;
plot(WN,V, 'k-', 'linewidth', 1);    % plot the W-nullcline
hold on;
plot(VN,V, 'k-.', 'linewidth', 1);   % plot the V-nullcline
hold on 
plot(Wss,Vss, '*', 'color','g');
hold on 

options = odeset('RelTol',1e-6,'AbsTol',[1e-8 1e-8]);
[T,Y] = ode15s(f,(0:0.1:450),[-1.,1.5],options);
[T1,Y1] = ode15s(f1,(0:0.1:450),[0,0.5],options);
%[T Y] = ode45s(f, [0 150], [1, 0]);
plot(Y(:,2), Y(:,1), 'b-', 'linewidth', 1.5);
hold on
plot(Y1(:,2), Y1(:,1), 'r-', 'linewidth', 1.5);
title(['Phase plane, I=', num2str(I)], 'fontsize', 16);
xlabel('V(t)'); ylabel('W(t)');
axis([-.8 2. -2.5 2.5]);

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  2     %%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2); hold off;
plot(T, Y(:,1), 'b-', 'linewidth', 1.5); hold on;
plot([0 90],[Vss Vss],'k:', 'linewidth', 1.5);
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
title(['Real(\lambda)=',num2str(real(Lambda(1)))], 'fontsize', 16);
axis tight;

%%%%%%%%%%%%%%%%%%%%%%%     PANEL  4     %%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4); hold on;
[T Y] = ode15s(f, [0 400], Y(length(T),:) );  % Run some more in order
[T1 Y1] = ode15s(f1, [0 400], Y1(length(T),:) );  % to settle at equilibrium
plot(I, min(Y(:,1)), 'b.');       % minimum of V(t): a magenta point
plot(I, max(Y(:,1)), 'b.');       % maximum of V(t): a magenta point
plot(I,Vss,'color','b');    % Fixed point: a red or blue point
plot(I, min(Y1(:,1)), 'r.');       % minimum of V(t): a magenta point
plot(I, max(Y1(:,1)), 'r.');       % maximum of V(t): a magenta point
title('Bifurcation diagram (V_{min}, V_{max}, V*)', 'fontsize', 12);
xlabel('I'); ylabel('V equilibrium');
axis tight;
drawnow;

end;



