function OrbitPlot(Y,Ye,RV_Earth,RV_Ast,X0min,~,trajmin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotates spacecraft state at Aphelion into the Sun-Centred inertial frame 
% and plots it's trajectory and asteroid trajectory about Sun for 2 years.
%                                                         
% INPUTS:
% - System of first order ODEs, dXdt
%
% OUTPUTS:                                               
% - Integrated velocity components, Position (km)
% - Integrated acceleration components, Velocity (km/s)
% - Plot of position and velocity Absolute Percentage Errors (APE) (%)
% - Plot of the norm of specific angular momentum, norm(h) (km^2/s)
%
% Author: Sho Wright (sw01745)                           
% Date created: 25/10/22                                         
% Date updated: 18/11/22                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

choice  = input(['Choose from the following: \n' ...
                 '1 for Earth and S/C escape from SEL2 to Aphelion: \n' ...
                 '2 for Sun and S/C transfer trajectory to asteroid:\n' ...
                 '3 for both diagrams:']);

%% Rotating Spacecraft into Inertial Frame
AU = 1.4959787070e08; mu = 1.32712440018e11; n = sqrt(mu/AU^3); 
mubar = 3.0032080443e-06;
arr = trajmin(1,:); dep = trajmin(2,:);
dt = trajmin(6,:);


r_sci = zeros(length(RV_Earth),3);
v_sci = zeros(length(RV_Earth),3);
for i = 1:length(RV_Earth)
    [r_sci(i,:),v_sci(i,:)] = Rotation(Ye(end,:),RV_Earth(i,:),mubar); %Rotates starting position at Perihelion into inertial frame
end
r_sci = r_sci*AU; v_sci = v_sci*(AU*n);
Xsc = [r_sci v_sci];

%Escape trajectory from SEL2 to Aphelion
r_esc = Y(:,1:3)*AU; v_esc = Y(:,4:6)*(AU*n);
X_esc = [r_esc v_esc];

Z = RV_Ast; %Asteroid State

% % Returning Orbital Elements for Initial States
% fprintf('Spacecraft at Perhelion entry \n')
% SCcoe= RV2COE(r_sci(i,:)', v_sci(i,:)',mu); 
% save('C:\Users\sho_w\Documents\University of Surrey\EEEM004 - Dissertation\Scripts\PorkChop\SCcoe.mat','SCcoe')

switch choice 
    case 1
        %% Plotting Earth and S/C Transfer
        R_Earth = 6378e3;
        Msun = 1.989e30; Mearth =5.972e24; AU = 1.4959787070e08;
        R_SOI = 0.9431*AU*(Mearth/Msun)^(2/5);
        %Pulling in surface wrap
        EarthWrap = imread('EarthWrap.jpg');
        [Xe,Ye,Ze] = sphere(50);
        % Create figure
        figure('units','normalized','outerposition',[0 0 1 1])
        % Plot Earth
        globe = surf(R_Earth*Xe, R_Earth*Ye, R_Earth*Ze, 'EdgeColor', 'none', 'FaceColor', 'none'); hold on; axis equal;
        globe2 = surf(R_SOI*Xe, R_SOI*Ye, R_SOI*Ze, 'EdgeColor', 'r', 'FaceColor', 'none'); hold on; axis equal;
        set(globe, 'FaceColor', 'texturemap', 'CData', EarthWrap, 'FaceAlpha', 0.9, 'EdgeColor', 'none');
        set(globe2, 'FaceColor', 'texturemap','FaceAlpha', 0.9, 'EdgeColor', 'none');
        xlabel('$X$, SCI (km)');
        ylabel('$Y$, SCI (km)');
        zlabel('$Z$, SCI (km)');
        title('SEL2-Aphelion: Transfer Design 1')
        % Plot ECI axes
        quiver3(0, 0, 0, 1, 0, 0, 1e4, 'k', 'LineWidth', 2);    % I-axis
        quiver3(0, 0, 0, 0, 1, 0, 1e4, 'k', 'LineWidth', 2);    % J-axis
        quiver3(0, 0, 0, 0, 0, 1, 1e4, 'k', 'LineWidth', 2);    % K-axis
        text(1e4, 0, 0, 'I', 'FontSize', 15, 'Interpreter', 'tex', 'Color', 'k')
        text(0, 1e4, 0, 'J', 'FontSize', 15, 'Interpreter', 'tex', 'Color', 'k')
        text(0, 0, 1e4, 'K', 'FontSize', 15, 'Interpreter', 'tex', 'Color', 'k')
        % Plot Transfer Trajectory for spacecraft 
        plot3(X_esc(:,1), X_esc(:,2), X_esc(:,3),'--','Color', 'k','LineWidth', 2);
        plot3(X_esc(1,1), X_esc(1,2), X_esc(1,3), 'ok', 'MarkerFaceColor', 'b');
        plot3(X_esc(end,1), X_esc(end,2), X_esc(end,3), 'ok', 'MarkerFaceColor', 'r');
        % Plot Trajectory for spacecraft
        plot3(Xsc(200:350,1), Xsc(200:350,2), Xsc(200:350,3),'k','LineWidth', 1.5);
        legend('','','','SCI Frame','Transfer Trajectory','SEL2','Heliocentric Entry','Spacecraft Orbit','Asteroid Orbit')
        az = -45;
        el = 45;
        view(az, el);

    case 2
        %% ODE Function for Transfer
        X0 = [X0min(1:3,:)*AU ; X0min(4:6,:)*AU*n];
        %Xf = [Xfmin(1:3,:)*AU ; Xfmin(4:6,:)*AU*n];
        t0 = 0; tspan = [t0 dt];
        options = odeset('RelTol',3e-14,'AbsTol',1e-16); %Setting tolerances (varied for situation)
        [~,T] = ode45(@(t,X)EOM2BP(t,X,mu),tspan,X0',options); 

        %% Plotting x-z View
        figure; 
        plot(Xsc(:,2),Xsc(:,3),'LineWidth', 2)
        hold on
        plot(Z(:,2),Z(:,3),'LineWidth', 2);
        title ('X-Z axis view of Non-Co-Planar Trajectories')
        xlabel ('SCI x-axis (km)')
        xlim([-2e08 2e08])
        ylabel ('SCI z-axis (km)')
        ylim([-8e07 8e07])
        box on; grid on; 
        legend('Spacecraft Orbit','Asteroid Orbit')

        %% Plotting Sun and Earth 
        R_Sun = 696000;%*1e01; %Radius of Sun (km)
        R_Earth = 6378;%*5e02;
        %Pulling in surface wrap
        SunWrap = imread('SunWrap.jpg');
        EarthWrap = imread('EarthWrap.jpg');
        [Xe,Ye,Ze] = sphere(50);
        % Create figure
        figure('units','normalized','outerposition',[0 0 1 1])
        % Plot Earth
        globe = surf(R_Sun*Xe, R_Sun*Ye, R_Sun*Ze, 'EdgeColor', 'none', 'FaceColor', 'none'); hold on; axis equal;
        globe2 = surf(R_Earth*Xe + AU, R_Earth*Ye, R_Earth*Ze, 'EdgeColor', 'none', 'FaceColor', 'none'); hold on; axis equal;
        set(globe2, 'FaceColor', 'texturemap', 'CData', EarthWrap, 'FaceAlpha', 0.9, 'EdgeColor', 'none');
        set(globe, 'FaceColor', 'texturemap', 'CData', SunWrap, 'FaceAlpha', 0.9, 'EdgeColor', 'none');
        xlabel('$X$, SCI (km)');
        ylabel('$Y$, SCI (km)');
        zlabel('$Z$, SCI (km)');
       
        % Plot ECI axes
        quiver3(0, 0, 0, 1, 0, 0, 0.5e8, 'k', 'LineWidth', 2);    % I-axis
        quiver3(0, 0, 0, 0, 1, 0, 0.5e8, 'k', 'LineWidth', 2);    % J-axis
        quiver3(0, 0, 0, 0, 0, 1, 0.5e8, 'k', 'LineWidth', 2);    % K-axis
        text(0.55e8, 0, 0, 'I', 'FontSize', 15, 'Interpreter', 'tex', 'Color', 'k')
        text(0, 0.55e8, 0, 'J', 'FontSize', 15, 'Interpreter', 'tex', 'Color', 'k')
        text(0, 0, 0.55e8, 'K', 'FontSize', 15, 'Interpreter', 'tex', 'Color', 'k')
        
        % Plot trajectory for spacecraft
        plot3(Xsc(1:dep,1), Xsc(1:dep,2), Xsc(1:dep,3),'--','Color','k','LineWidth', 1.5);
        % plot3(Xsc(1,1), Xsc(1,2), Xsc(1,3), 'ok', 'MarkerFaceColor', 'b');
        % plot3(Xsc(end,1), Xsc(end,2), Xsc(end,3), 'ok', 'MarkerFaceColor', 'r');
        hold on
        
        % Plot trajectory for asteroid initial
        plot3(Z(1:arr,1), Z(1:arr,2), Z(1:arr,3),':','Color', 'k','LineWidth', 1.5);
        % plot3(Z(1,1), Z(1,2), Z(1,3), 'ok', 'MarkerFaceColor', 'g');
        % plot3(Z(end,1), Z(end,2), Z(end,3), 'ok', 'MarkerFaceColor', 'm');

        % Plot transfer trajectory from spacecraft to asteroid
        plot3(T(:,1), T(:,2), T(:,3), 'r','LineWidth', 2);
        plot3(T(1,1), T(1,2), T(1,3), 'ok', 'MarkerFaceColor', 'b');
        plot3(T(end,1), T(end,2), T(end,3), 'ok', 'MarkerFaceColor', 'g');
        
        title('Transfer Trajectory to Asteroid (Minimum departure velocity)','FontSize',9)
        legend('','','','','SCI Frame','Spacecraft Orbit','Asteroid Orbit','Transfer Trajectory','Departure','Arrival','FontSize',9)
        az = -45;
        el = 45;
        view(az, el);
end
