function [L2, V] = Lagrange(mubar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using Newton.m script, sets initial guesses for each of the 3 collinear
% equilibrium points and plots them alongside the 2 triangular equilibrium
% points. 
%
% To test this function for the Sun-Earth System copy and paste:
%   - Lagrange(3e-06)
%
% INPUTS:
% - Mass-ratio parameter, mubar, (Unitless)
%   - mubar = m2/m1+m2
%   - where m1 and m2 are the masses of the primaries (kg)
%   - E.g., for a Sun-Earth system:
%        Mass of Sun:   m1 = 1.989E30 (kg)
%        Mass of Earth: m2 = 5.972E24 (kg)                                                    
% Other common mubar values:
% - Earth-Moon: 1.2151E-2
% - Sun-Jupiter: 7.1904E-4
% - Sun-Saturn: 2.8571E-4
% - Saturn-Titan: 2.366E-4
%
% OUTPUTS:                                               
% - 3 Collinear equilibrium points L1, L2, and L3
% - Final EOM values evaluated at L1, L2, and L3
% - 6D Eigenvector, stability properties of L2 (Unitless)
% - Plot of the 5 equilibrium points alongside the primaries
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 02/04/23                                         
% Date updated: 20/04/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting up script
close all; clc; format short;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 16);

%% Newton's method for L1
for X0 = [(1-mubar)*0.98;0;0;0;0;0] %Initial guess for L1
    [Xstar,i] = Newton(mubar,X0); %Applying Newton's method
    L1 = Xstar; %Defining final converged vector as L1
end

%% Newton's method for L2
for X0 = [(1-mubar)*1.02;0;0;0;0;0] %Initial guess for L2
    [Xstar,i] = Newton(mubar,X0); %Applying Newton's method
    L2 = Xstar; %Defining final converged vector as L2
end

%% Newton's method for L3
for X0 = [-1.2;0;0;0;0;0] %Initial guess for L3
    [Xstar,i] = Newton(mubar,X0); %Applying Newton's method
    L3 = Xstar; %Defining final converged vector as L3
end

%% Outputting iteration to convergence and X vector at convergence
fprintf('X converges to L1 after %.0f iterations \n',i) 
fprintf('X vector at L1:\n')
disp(transpose(L1))

fprintf('X converges to L2 after %.0f iterations \n',i)
fprintf('X vector at L2:\n')
disp(transpose(L2))

fprintf('X converges to L3 after %.0f iterations \n',i)
fprintf('X vector at L3:\n')
disp(transpose(L3))

%% Outputting L2 Eigenvector
J = Jacobi(mubar,L2); %Jacobian matrix evaluated at L2 
[~,V] = eig(J); %Returining 6D eigenvector of Jacobian matrix at L2

fprintf('Eigenvalues of final Jacobian matrix for L2: \n')
disp(eig(J))

%% Temporarily putting perturbation here for a figure
ExpGrowth = V(:,1); %Slicing out positive real component (exponential growth)
NonD = vecnorm(ExpGrowth);
ExpGrowth(:,1) = ExpGrowth(:,1)/NonD; %Normalising eigenvector
deltaX = ExpGrowth*3e-05; %Scaling factor 1e-06
X0 = L2 + deltaX; %New initial condition
X1 = L2 - deltaX;

%% Plotting figure
figure;
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); %Enlarge figure to full screen
title('SEL2 Escape for Design 2')
hold on
xlim([-L2(1,1)*1.1,L2(1,1)*1.1])
ylim([-L2(1,1)*1.1,L2(1,1)*1.1])
grid on 
axis equal 
xlabel('x (rotating frame)')
ylabel('y (rotating frame)')
plot(-mubar,0,'bo','MarkerFaceColor','b') %Plotting mass 1
plot(1-mubar,0,'go','MarkerFaceColor','g') %Plotting mass 2
plot(L1(1,1),0,'ksquare','MarkerFaceColor','k') %Plotting L1
plot(L2(1,1),0,'kpentagram','MarkerFaceColor','k') %Plotting L2
plot(L3(1,1),0,'khexagram','MarkerFaceColor','k') %Plotting L3

plot(X0(1,1),0,'rpentagram','MarkerFaceColor','r') %Plotting perturbation
plot(X1(1,1),0,'bpentagram','MarkerFaceColor','b') %Plotting perturbation

plot(0.5-mubar,sqrt(3)/2,'r^','MarkerFaceColor','r') %Plotting L4
plot(0.5-mubar,-sqrt(3)/2,'rv','MarkerFaceColor','r') %Plotting L5
plot([-mubar,0.5-mubar,1-mubar,0.5-mubar,-mubar], [0,sqrt(3)/2,0,-sqrt(3)/2,0], '--k')
yline(0,'k')
coords = linspace(0,pi,100); %Slicing range 0 to pi into 100 points
plot((1-mubar).*cos(coords),(1-mubar).*sin(coords),'Color',"#7E2F8E") %Plotting upper semicircle
hold all
plot((1-mubar).*cos(coords),-(1-mubar).*sin(coords),'Color',"#7E2F8E") %Plotting lower semicircle
legend('','', '', '$L_2$','','$IC_1$ (Design 1)','$IC_2$ (Design 2)','Interpreter','latex','FontSize',12)
% legend('$m_1$','$m_2$', '$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$','Interpreter','latex','FontSize',12)

%% Clearing 'ans' in Command Window
if nargout<1 
    clear L2 V 
end