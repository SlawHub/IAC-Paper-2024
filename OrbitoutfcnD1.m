function stop = OrbitoutfcnD1(Xopt,~,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots first mission trajectories (prior to optimisation) and the output
% trajectory of the last iteration of fmincon (optimised trajectory)
%
% INPUTS:
% - Optimisation vector X (14 x 1), including:
%   - State space at SEL2, X0 (6 x 1)
%   - Time of SEL2 escape, TOF1 (1 x 1)
%   - Departure velocity, DV1 (3 x 1)
%   - Time of transfer, TOF2 (1 x 1)
%   - Relative velocity, DV2 (3 x 1)
% - fmincon structure, optimvalues
% - fmincon iteration state flag, including:
%   - init: initialisation
%   - iter: all iterations preceding
%   - done: last iteration state
%   - interrupt: for errors
%
% OUTPUTS:                                               
% - Orbit plot diagram: Pre and Post optimised trajectories
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 12/08/23                                         
% Date updated: 12/08/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining Variables
mubar = 3.0032080443e-06; %Mass-ratio parameter of Sun-Earth
muSun = 1.32712440018e11; 
LU = 1.4959787070e08; %Length Unit
n = sqrt(muSun/LU^3); %Mean Motion of Sun
TU = 1/n; %Time Unit
JDJ2000 = 2451545.0;
sidereal = 86400;
%Pulling in planet wraps
R_Sun = 696340; %Radius of Sun (km)
SunWrap = imread('SunWrap.jpg');
mu = muSun*TU^2/LU^3; 
stop = false; %This is just needed, don't ask me why

%% Plotting Orbit
switch flag
    case 'init'    
    case 'iter'
 
    case 'done'
        %% Loading in D1 CR3BP
        load('D1CR3BP.mat','t1e');
        TOF1ori = t1e*TU;

        %% Retrieving Asteroid State
        rmpath(genpath('AsteroidSearchResults'))
        Design = 1; spkid = 3654379; 
        DesignFilePath = sprintf('AsteroidSearchResults/D%d',Design);
        addpath(genpath(DesignFilePath))
        bspfilename = sprintf('%d.bsp',spkid); 
        cspice_furnsh(which(bspfilename));
        datei = datetime(2026,1,1,00,00,00);
        datef = datetime(2029,5,15,00,00,00);
        DATEarr = datei:days(1):datef;
        JD_Arr = juliandate(DATEarr);
        ephemeris_time = JD_Arr*sidereal - JDJ2000*sidereal;
        RV_Ast = cspice_spkezr(sprintf('%d',spkid),ephemeris_time,'ECLIPJ2000','NONE','10');  

        %% Optimised escape trajectory from SEL2 to Aphelion
        options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16); %Setting tolerances (varied for situation)
        [TOF1opt,Y1] = ode45(@(t,X)CR3BPEOM(mubar,X,t),[0 Xopt(7)],Xopt(1:6),options); %Integrating system of 1st order ODEs
        
        %% Finding Real SEL2 Departure Date
        NEODepDate = datetime(Xopt(15),'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');
        NEOetd = cspice_str2et(cellstr(NEODepDate));
        SEL2DepDate = NEODepDate - duration(seconds(Xopt(7)*TU));
        SEL2etd = cspice_str2et(cellstr(SEL2DepDate));
        etds = linspace(SEL2etd,NEOetd,length(Y1));
        
        %% Rotating Y1 trajectory into ECLJ2K frame
        XLeg1 = zeros(6,length(Y1));
        for i = 1:1:length(etds)
            X1endSCS = (Y1(i,:)' + [mubar 0 0 0 0 0]'); %Aphelion in normalised Sun-Centered Synodic frame 
            X1endSCS = [X1endSCS(1:3)*LU ; X1endSCS(4:6)*LU/TU]; %Aphelion in SCS frame dimensionalised
            STM = cspice_sxform('SUNEARTH','ECLIPJ2000', etds(i)); %State Transformation Matrix
            XLeg1(:,i) = STM*X1endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        end
        XLeg1 = XLeg1';

        %% Adding Burn
        %Using Spice syntax
        Time = Xopt(15)*sidereal - TOF1ori + TOF1opt(end)*TU;
        Time = Time - 2451545.0*sidereal; %2451545.0 = JD 1st Jan 
        %Rotating end of Leg 1 (Aphelion) into SCI frame
        Y1endSCS = (Y1(end,:)' + [mubar 0 0 0 0 0]'); %Aphelion in normalised Sun-Centered Synodic frame 
        Y1endSCS = [Y1endSCS(1:3)*LU ; Y1endSCS(4:6)*LU/TU];
        STM = cspice_sxform('SUNEARTH','ECLIPJ2000', Time); %State Transformation Matrix  
        Y1 = STM*Y1endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        Y1new = [Y1(1:3)/LU ; Y1(4:6)/LU*TU];
     
        %% Propagating Leg 2 
        options = odeset('RelTol',3e-14,'AbsTol',1e-16); %Setting tolerances (varied for situation)
        [~,Y2] = ode45(@(t,X)EOM2BP(t,X,mu),[0 Xopt(11)],Y1new',options); 
        XLeg2 = [Y2(:,1:3)*LU Y2(:,4:6)*LU/TU];

% %         %Rotating end of Leg 1 (Aphelion) into SCI frame
%         YLeg1end = [Leg1(end,1:3)/LU , Leg1(end,4:6)/LU*TU];
%         Y1endSCS = (YLeg1end' + [mubar 0 0 0 0 0]'); %Aphelion in normalised Sun-Centered Synodic frame 
%         Y1endSCS = [Y1endSCS(1:3)*LU ; Y1endSCS(4:6)*LU/TU];
%         STM = cspice_sxform('SUNEARTH','ECLIPJ2000', Time); %State Transformation Matrix  
%         Y2 = STM*Y1endSCS; %Rotating into SCI (Ecliptic J2000 frame)
%         Y2 = [Y2(1:3)/LU ; Y2(4:6)/LU*TU];
%         Y2new = Y2 + [zeros(3,1) ; Xopt(8:10)]; %Transfer state space vector

        %% Plotting Full Initially Found Mission
        % Create figure
        fig1 = figure(1);
        % Plot Sun
        [Xe,Ye,Ze] = sphere(50);
        globe = surf(R_Sun*Xe, R_Sun*Ye, R_Sun*Ze, 'EdgeColor', 'none', 'FaceColor', 'none'); 
        hold on; axis equal;
        set(globe, 'FaceColor', 'texturemap', 'CData', SunWrap, 'FaceAlpha', 0.9, 'EdgeColor', 'none');
        xlabel('$X$, ECLJ2K (km)');
        ylabel('$Y$, ECLJ2K (km)');
        zlabel('$Z$, ECLJ2K (km)');  
        % Plot ECI axes
        quiver3(0, 0, 0, 1, 0, 0, 0.5e8, 'k', 'LineWidth', 2);    % I-axis
        quiver3(0, 0, 0, 0, 1, 0, 0.5e8, 'k', 'LineWidth', 2);    % J-axis
        quiver3(0, 0, 0, 0, 0, 1, 0.5e8, 'k', 'LineWidth', 2);    % K-axis
        text(0.55e8, 0, 0, 'I', 'FontSize', 15, 'Interpreter', 'tex', 'Color', 'k')
        text(0, 0.55e8, 0, 'J', 'FontSize', 15, 'Interpreter', 'tex', 'Color', 'k')
        text(0, 0, 0.55e8, 'K', 'FontSize', 15, 'Interpreter', 'tex', 'Color', 'k') 
        
        %Plot SEL2 escape Leg 1 
        plot3(XLeg1(1,1), XLeg1(1,2), XLeg1(1,3), 'ok', 'MarkerFaceColor', 'b','MarkerSize',8); %4000:end
        plot3(XLeg1(:,1), XLeg1(:,2), XLeg1(:,3),'--k','LineWidth', 3.5);
        plot3(XLeg1(end,1), XLeg1(end,2), XLeg1(end,3), 'ok', 'MarkerFaceColor', 'g','MarkerSize',8);
        
        %Plot transfer for Leg 2
        % plot3(XLeg2(1,1), XLeg2(1,2), XLeg2(1,3), 'ok', 'MarkerFaceColor', 'g','MarkerSize',8);
        plot3(XLeg2(:,1), XLeg2(:,2), XLeg2(:,3),'-k','LineWidth', 3.5);
        plot3(XLeg2(end,1), XLeg2(end,2), XLeg2(end,3), 'ok', 'MarkerFaceColor', 'r','MarkerSize',8);
        
        %Plot orbit for asteroid 
        plot3(RV_Ast(1,:), RV_Ast(2,:), RV_Ast(3,:),'-','Color', 'b','LineWidth', 3);
        
        axis equal 
        ax = gca;
        ax.FontSize = 14; 
        ax.FontSize = 14; 
        
        title('Mission Scenario 1','FontSize',18)
%         legend('','','','SCI Frame','SEL2','Manifold Transfer','Aphelion','NEO Transfer','Rendezvous','NEO Orbit','FontSize',16);
%         exportgraphics(fig1,'D12013XV8FINALFULLLEGEND.png','Resolution',1000)
end



% % Set last_best to best
% if optimvalues.iteration == 0
% last_best = best;
% else
%     plot3(Y2opt(:,1), Y2opt(:,2), Y2opt(:,3),':','Color', 'b','LineWidth', 1.5);
% end
