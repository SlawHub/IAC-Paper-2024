function [Ye,te] = CR3BPMain(mubar,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By perturbing the spacecraft by the positive real eigenvalue (the 
% exponential growth component), it is possible to 'nudge' the s/c into 
% the exterior realm (deep space) and use an ode function to propagate it 
% along the unstable manifold. From here an event function can be applied 
% to stop integration once the s/c reaches heliocentric perihelion. 
%
% - To test: CR3BPMain(3.0032080443e-06,50)
%
% INPUTS:
% - Mass-ratio parameter, mubar, (Unitless)			                                                     
% - Propagation time, dt (hours)
% 
% OUTPUTS:                                               
% - Final position and velocity components of satellite in rotating
% reference frame at entry into Heliocentric orbit at the Perihelion, Ye
% - Time of event, te
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 21/03/23                                         
% Date updated: 19/04/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting up script
close all; clc; format short;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 16);

%% Defining Variables
sidereal = 86400; AU = 1.4959787070e08; 
Smu = 1.32712440018e11; Emu = 3.986004418e05;
n = sqrt((Smu + Emu)/AU^3); 

%% User Input
choice = input('Choose 1 for for design 1 (via unstable manifold), or 2 for design 2 (via powered flyby):');

%% Determining new intial condition 
[L2, V] = Lagrange(mubar); %Pulling in X for Collinear Lagrange Points
ExpGrowth = V(:,1); %Slicing out positive real component (exponential growth)
NonD = vecnorm(ExpGrowth);
ExpGrowth(:,1) = ExpGrowth(:,1)/NonD; %Normalising eigenvector
deltaX = ExpGrowth*3e-05; %Scaling factor 1e-06 1e-12

switch choice
    case 1
        X0 = L2 + deltaX; %New initial condition
        options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16,'OutputFcn',@odephas2,'Events',@HelioEvent); %Setting tolerances (varied for situation)

    case 2
        X0 = L2 - deltaX;
        options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16,'OutputFcn',@odephas2,'Events',@GeoEvent); %Setting tolerances (varied for situation)

end

%% Integrating CR3BP EOM
tspan = [0 dt]; %Defining time vector
[t,Y,te,Ye] = ode45(@(t,X)EOM(mubar,X,t),tspan,X0',options); %Integrating system of 1st order ODEs
%Adding 3 new outputs: time of event, vector of event, iteration of event

%% Creating new vector adding 150m/s impulsive burn
switch choice
    case 2
%         Emu = Emu/AU^3/n^2;
        fuel = 150e-03/AU/n;
        vnew = Ye(end,4:6) + fuel;
        X0 = [Ye(end,1:3), vnew ];
        options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16,'OutputFcn',@odephas2,'Events',@HelioEvent); %Setting tolerances (varied for situation)
        [t2,Z,te2,Ze] = ode45(@(t,X)EOM(mubar,X,t),tspan,X0',options); %Integrating system of 1st order ODEs
%         options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16,'OutputFcn',@odephas2,'Events',@SOIEvent); %Setting tolerances (varied for situation)        
%         [t2,Z,te2,Ze] = ode45(@(t,X)EOM2BP(t,X,Emu),tspan,X0',options); %Integrating system of 1st order ODEs

end 

%Starting code for velocity globe:
% lon = linspace(0,2*pi,100);
% lat = linspace(-pi,pi,100);
% vglobe = meshgrid(lon,lat);

%% Plotting Output
figure; 
switch choice
    case 1
        plot(Y(:,1),Y(:,2),Ye(:,1),Ye(:,2),'o','LineWidth',1.5);
        hold on
        title ('SEL2 to Aphelion - via SEL2 Homoclinic Orbit')
    case 2 
        plot(Y(:,1),Y(:,2),Ye(:,1),Ye(:,2),'o','LineWidth',1.5);
        hold on
        plot(Z(:,1),Z(:,2),Ze(:,1),Ze(:,2),'o','LineWidth',1.5);
        plot(1-mubar,0,'ko','MarkerFaceColor','r') %Plotting mass 2
        title('SEL2 to Aphelion - via Earth Powered Flyby') 
        legend('SEL2 to Perigee','','Perigee to Aphelion','FontSize',10)
end
hold off
ylabel ('y(t) (rotating frame)')
xlabel ('x(t) (rotating frame)')
box on; grid on;

%% Printing Results

switch choice
    case 1
        Tescape = te(end,:)*n^-1;
        fprintf('Time to escape SEL2 to the Aphelion: %.3f days \n',Tescape/sidereal);
        fprintf('Vector at entry into Aphelion in rotating frame:\n')
        disp(Ye(end,:))
    case 2
        Tescape1 = te(end,:)*n^-1;
        Tescape2 = te2(end,:)*n^-1;
        fprintf('Time to escape SEL2 to the Perigee: %.3f days \n',Tescape1/sidereal);
        fprintf('Time to transfer from Perigee to Aphelion: %.3f days \n',Tescape2/sidereal);
        fprintf('Total Earth powered flyby maneouvre time: %.3f days \n', (Tescape1+Tescape2)/sidereal);
        fprintf('Vector at entry into Perigee in rotating frame:\n')
        disp(Ye(end,:))
        fprintf('Vector at entry into Aphelion in rotating frame:\n')
        disp(Ze(end,:))
end

%% Calculating and Plotting Jacobi Constant
figure;
switch choice 
    case 1
        C = Constant(mubar,Y); %Calling Jacobi constant script
        plot(t,C,'r','linewidth',0.1);
        xlim([0,te(end,:)]); ylim([1.0, 1.25*C(end,:)]);
        legend('C');
    case 2
        C1 = Constant(mubar,Y); %Calling Jacobi constant script
        C2 = Constant(mubar,Z); %Calling Jacobi constant script
        plot(t,C1,'r','linewidth',0.1);
        hold on
        plot(t2,C2,'b','linewidth',0.1);
        xlim([0,te(end,:)]); ylim([1.0, 1.25*C1(end,:)]);
        legend('C1','C2');
end
title('Jacobi Constant throughout Transfer');
ylabel('Jacobi Constant (-)','Interpreter','none');
xlabel('Time');
box on; grid on;

%% Storing Outputs
switch choice
    case 1
        save('C:\Users\sho_w\Documents\University of Surrey\EEEM004 - Dissertation\Scripts\PorkChop\CR3BPOutput1.mat','Y','Ye','te')
        save('C:\Users\sho_w\Documents\University of Surrey\EEEM004 - Dissertation\Scripts\MOID\CR3BPOutput1.mat','Y','Ye','te')
    case 2
%         save('C:\Users\sho_w\Documents\University of Surrey\EEEM004 - Dissertation\Scripts\PorkChop\CR3BPOutput2.5.mat','Y','Ye','te')
%         save('C:\Users\sho_w\Documents\University of Surrey\EEEM004 - Dissertation\Scripts\MOID\CR3BPOutput2.5.mat','Y','Ye','te')  
        Y = Z; Ye = Ze;
        te = te(end,:) + te2(end,:);
        save('C:\Users\sho_w\Documents\University of Surrey\EEEM004 - Dissertation\Scripts\PorkChop\CR3BPOutput2.mat','Y','Ye','te')
        save('C:\Users\sho_w\Documents\University of Surrey\EEEM004 - Dissertation\Scripts\MOID\CR3BPOutput2.mat','Y','Ye','te')
end
%% Clearing 'ans' in Command Window
if nargout<1 
    clear Ye te
end

%% Plotting Homoclinic Orbits
% X0 = L2 + deltaX;
% X1 = L2 - deltaX;
% options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16,'OutputFcn',@odephas2,'Events',@HelioEvent); %Setting tolerances (varied for situation)
% dt = 95.6;
% tspan = [0 dt];
% [~,Y,~,~] = ode45(@(t,X)EOM(mubar,X,t),tspan,X0',options); %Integrating system of 1st order ODEs
% dt = 139.8;
% tspan = [0 dt];
% [~,Z,~,~] = ode45(@(t,X)EOM(mubar,X,t),tspan,X1',options); %Integrating system of 1st order ODEs
% 
% figure;
% % plot(Y(:,1),Y(:,2),'Color',"#0072BD",'LineWidth',1.5);
% hold on
% % plot(Z(:,1),Z(:,2),'Color',"#77AC30",'LineWidth',1.5);
% plot(Z(:,1),Z(:,2),'Color',"#4DBEEE",'LineWidth',1.5);
% coords = linspace(0,pi,100); %Slicing range 0 to pi into 100 points
% plot((1-mubar).*cos(coords),(1-mubar).*sin(coords),'--','Color',"#7E2F8E",'LineWidth',1.5) %Plotting upper semicircle
% hold all
% plot((1-mubar).*cos(coords),-(1-mubar).*sin(coords),'--','Color',"#7E2F8E",'LineWidth',1.5) %Plotting lower semicircle
% p = plot(-mubar,0,'ko','MarkerFaceColor','#D95319'); %Plotting mass 1
% p.MarkerSize = 8;
% q = plot(1-mubar,0,'ko','MarkerFaceColor','r'); %Plotting mass 2
% q.MarkerSize = 4;
% % title('Homoclinic Orbits')
% title('Heteroclinic Connection from SEL2 to SEL1')
% hold off
% ylabel ('y(t) (rotating frame)')
% xlabel ('x(t) (rotating frame)')
% box on; grid on;
% % legend('SEL2','SEL1','FontSize',10)