function stop = OrbitoutfcnD2(Xopt,~,flag)
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
        %% Retrieving Earth Data
        %Earth state in day resolution
        fid = fopen('Earth.txt'); 
        A = textscan(fid,'%f %s %f %f %f %f %f %f', 'Delimiter', ',')';
        fclose(fid);
        DATEdep = A{2}';
        %Slicing out tick labels for each axis
        DATEdep = cellfun(@(x)x(1:end-13),DATEdep,'un',0);
        DATEdep = cellfun(@(x)x(6:end),DATEdep,'un',0);
        etd = cspice_str2et(DATEdep); %Ephemeris time at Aphelion

        %% Loading in D1 CR3BP
        load('D2CR3BP.mat','t1e');
        TOF1ori = t1e*TU;
        load('D:\Optimisation\AsteroidSearchResults\D2\Transfer_Files\D2Transfer469219 Kamo`oalewa (2016 HO3).mat','trajmin');
        etPerigee = trajmin(8);
        dep = trajmin(2);
        TOF2ori = (etd(dep) - etPerigee)/TU;

        %% Retrieving Asteroid State
        rmpath(genpath('AsteroidSearchResults'))
        Design = 2; spkid = 54298612;                                                 % CHANGE ME 
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

        %% Generating Full Mission 
        %Leg 1 - Escaping SEL2 (All designs)
        options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16); %Setting tolerances (varied for situation)
        [TOF1opt,Y1opt] = ode45(@(t,X)CR3BPEOM(mubar,Xopt,t),[0 Xopt(7)],Xopt(1:6),options); %Integrating system of 1st order ODEs
        TOF1opt = TOF1opt(end)*TU; %Dimensionalise
        %Rotating DV1 to Ecliptic J2K frame and adding for EPF
        DV1EJ2K = [zeros(3,1)*LU; Xopt(8:10)*LU/TU];
        etd1 = Xopt(19)*sidereal - TOF1ori + TOF1opt;% - 2451545.0*sidereal;
        DateDep = datetime(etd1/sidereal,'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');
        etd1 = cspice_str2et(cellstr(DateDep));
        STM = cspice_sxform('ECLIPJ2000','SUNEARTH', etd1); %State Transformation Matrix 
        DV1Syn = DV1EJ2K'*STM; %DV1 in Synodic frame
        DV1Syn = [DV1Syn(1:3)/LU , DV1Syn(4:6)/LU*TU];
        Y1new = Y1opt(end,:)' + DV1Syn';
        %Propagating to SOI in CR3BP Domain
        [TOF2opt,Y2opt] = ode45(@(t,X)CR3BPEOM(mubar,X,t),[0 Xopt(11)],Y1new,options); %Integrating system of 1st order ODEs
        TOF2opt = TOF2opt(end)*TU; 
        %Rotating to Ecliptic J2000 Frame
        Time = Xopt(19)*sidereal - TOF1ori + TOF1opt - TOF2ori + TOF2opt;% - 2451545.0*sidereal;
        Time = Time - 2451545.0*sidereal;
        Y2endSCS = (Y2opt(end,:)' + [mubar 0 0 0 0 0]'); %SOI Edge in normalised Sun-Centered Synodic frame 
        Y2endSCS = [Y2endSCS(1:3)*LU ; Y2endSCS(4:6)*LU/TU];
        STM = cspice_sxform('SUNEARTH','ECLIPJ2000', Time); %State Transformation Matrix  
        Y2 = STM*Y2endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        Y2 = [Y2(1:3)/LU ; Y2(4:6)/LU*TU];
        % Final Tranfer to NEO
        X2 = Y2 + [zeros(3,1) ; Xopt(12:14)]; %Transfer state space vector
        [~, Y3opt] = ode45(@(t,X)EOM2BP(t,X,mu),[0 Xopt(15)],X2,options); %Integrating system of 1st order ODEs
        XLeg3 = [Y3opt(:,1:3)*LU Y3opt(:,4:6)*LU/TU];

        %% Getting Everything Into SCI Frame
        %Finding major dates
        NEOArrJD = Xopt(19)*sidereal + (Xopt(11) + Xopt(15))*TU;
        NEOArrDate = datetime(NEOArrJD/sidereal,'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');
        SOIDepDate = NEOArrDate - duration(seconds(Xopt(15)*TU));
        EPFDepDate = datetime(Xopt(19),'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');
        SEL2DepDate = EPFDepDate - duration(seconds(Xopt(7)*TU));
        %Rotating Y1 trajectory into ECLJ2K frame
        SEL2Depetd = cspice_str2et(cellstr(SEL2DepDate));
        EPFDepetd = cspice_str2et(cellstr(EPFDepDate));
        etds = linspace(SEL2Depetd,EPFDepetd,length(Y1opt));
        XLeg1 = zeros(6,length(Y1opt));
        for i = 1:1:length(etds)
            X1endSCS = (Y1opt(i,:)' + [mubar 0 0 0 0 0]'); %Aphelion in normalised Sun-Centered Synodic frame 
            X1endSCS = [X1endSCS(1:3)*LU ; X1endSCS(4:6)*LU/TU]; %Aphelion in SCS frame dimensionalised
            STM = cspice_sxform('SUNEARTH','ECLIPJ2000', etds(i)); %State Transformation Matrix
            XLeg1(:,i) = STM*X1endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        end
        XLeg1 = XLeg1';
        %Rotating Y2 trajectory into ECLJ2K frame
        SOIDepetd = cspice_str2et(cellstr(SOIDepDate));
        etds = linspace(EPFDepetd,SOIDepetd,length(Y2opt));
        %Rotating Y1 trajectory into ECLJ2K frame
        XLeg2 = zeros(6,length(Y2opt));
        for i = 1:1:length(etds)
            X1endSCS = (Y2opt(i,:)' + [mubar 0 0 0 0 0]'); %Aphelion in normalised Sun-Centered Synodic frame 
            X1endSCS = [X1endSCS(1:3)*LU ; X1endSCS(4:6)*LU/TU]; %Aphelion in SCS frame dimensionalised
            STM = cspice_sxform('SUNEARTH','ECLIPJ2000', etds(i)); %State Transformation Matrix
            XLeg2(:,i) = STM*X1endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        end
        XLeg2 = XLeg2';
        

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
        plot3(XLeg2(:,1), XLeg2(:,2), XLeg2(:,3),'-k','LineWidth', 3.5);
        plot3(XLeg2(end,1), XLeg2(end,2), XLeg2(end,3), 'ok', 'MarkerFaceColor', 'y','MarkerSize',8);
        %Plot transfer for Leg 3
        plot3(XLeg3(:,1), XLeg3(:,2), XLeg3(:,3),'-k','LineWidth', 3.5);
        plot3(XLeg3(end,1), XLeg3(end,2), XLeg3(end,3), 'ok', 'MarkerFaceColor', 'r','MarkerSize',8);
        %Plot orbit for asteroid 
        plot3(RV_Ast(1,:), RV_Ast(2,:), RV_Ast(3,:),'-','Color', 'b','LineWidth', 3);
        
        axis equal 
        ax = gca;
        ax.FontSize = 14; 
        ax.FontSize = 14; 
        
        title('Mission Scenario 1','FontSize',18)
        legend('','','','SCI Frame','SEL2','Manifold Transfer','Perigee','SOI Escape','SOI','NEO Transfer','Rendezvous','NEO Orbit','FontSize',16);
        exportgraphics(fig1,'D22016HO3FINALFULLLEGEND.png','Resolution',1000)
    case 'iter'
 
    case 'done'
        %% Retrieving Earth Data
        %Earth state in day resolution
        fid = fopen('Earth.txt'); 
        A = textscan(fid,'%f %s %f %f %f %f %f %f', 'Delimiter', ',')';
        fclose(fid);
        DATEdep = A{2}';
        %Slicing out tick labels for each axis
        DATEdep = cellfun(@(x)x(1:end-13),DATEdep,'un',0);
        DATEdep = cellfun(@(x)x(6:end),DATEdep,'un',0);
        etd = cspice_str2et(DATEdep); %Ephemeris time at Aphelion

        %% Loading in D1 CR3BP
        load('D2CR3BP.mat','t1e');
        TOF1ori = t1e*TU;
        load('D:\Optimisation\AsteroidSearchResults\D2\Transfer_Files\D2Transfer469219 Kamo`oalewa (2016 HO3).mat','trajmin');
        etPerigee = trajmin(8);
        dep = trajmin(2);
        TOF2ori = (etd(dep) - etPerigee)/TU;

        %% Retrieving Asteroid State
        rmpath(genpath('AsteroidSearchResults'))
        Design = 2; spkid = 54298612;                                                 % CHANGE ME 
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

        %% Generating Full Mission 
        %Leg 1 - Escaping SEL2 (All designs)
        options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16); %Setting tolerances (varied for situation)
        [TOF1opt,Y1opt] = ode45(@(t,X)CR3BPEOM(mubar,Xopt,t),[0 Xopt(7)],Xopt(1:6),options); %Integrating system of 1st order ODEs
        TOF1opt = TOF1opt(end)*TU; %Dimensionalise
        %Rotating DV1 to Ecliptic J2K frame and adding for EPF
        DV1EJ2K = [zeros(3,1)*LU; Xopt(8:10)*LU/TU];
        etd1 = Xopt(19)*sidereal - TOF1ori + TOF1opt;% - 2451545.0*sidereal;
        DateDep = datetime(etd1/sidereal,'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');
        etd1 = cspice_str2et(cellstr(DateDep));
        STM = cspice_sxform('ECLIPJ2000','SUNEARTH', etd1); %State Transformation Matrix 
        DV1Syn = DV1EJ2K'*STM; %DV1 in Synodic frame
        DV1Syn = [DV1Syn(1:3)/LU , DV1Syn(4:6)/LU*TU];
        Y1new = Y1opt(end,:)' + DV1Syn';
        %Propagating to SOI in CR3BP Domain
        [TOF2opt,Y2opt] = ode45(@(t,X)CR3BPEOM(mubar,X,t),[0 Xopt(11)],Y1new,options); %Integrating system of 1st order ODEs
        TOF2opt = TOF2opt(end)*TU; 
        %Rotating to Ecliptic J2000 Frame
        Time = Xopt(19)*sidereal - TOF1ori + TOF1opt - TOF2ori + TOF2opt;% - 2451545.0*sidereal;
        Time = Time - 2451545.0*sidereal;
        Y2endSCS = (Y2opt(end,:)' + [mubar 0 0 0 0 0]'); %SOI Edge in normalised Sun-Centered Synodic frame 
        Y2endSCS = [Y2endSCS(1:3)*LU ; Y2endSCS(4:6)*LU/TU];
        STM = cspice_sxform('SUNEARTH','ECLIPJ2000', Time); %State Transformation Matrix  
        Y2 = STM*Y2endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        Y2 = [Y2(1:3)/LU ; Y2(4:6)/LU*TU];
        % Final Tranfer to NEO
        X2 = Y2 + [zeros(3,1) ; Xopt(12:14)]; %Transfer state space vector
        [~, Y3opt] = ode45(@(t,X)EOM2BP(t,X,mu),[0 Xopt(15)],X2,options); %Integrating system of 1st order ODEs
        XLeg3 = [Y3opt(:,1:3)*LU Y3opt(:,4:6)*LU/TU];

        %% Getting Everything Into SCI Frame
        %Finding major dates
        NEOArrJD = Xopt(19)*sidereal + (Xopt(11) + Xopt(15))*TU;
        NEOArrDate = datetime(NEOArrJD/sidereal,'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');
        SOIDepDate = NEOArrDate - duration(seconds(Xopt(15)*TU));
        EPFDepDate = datetime(Xopt(19),'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');
        SEL2DepDate = EPFDepDate - duration(seconds(Xopt(7)*TU));
        %Rotating Y1 trajectory into ECLJ2K frame
        SEL2Depetd = cspice_str2et(cellstr(SEL2DepDate));
        EPFDepetd = cspice_str2et(cellstr(EPFDepDate));
        etds = linspace(SEL2Depetd,EPFDepetd,length(Y1opt));
        XLeg1 = zeros(6,length(Y1opt));
        for i = 1:1:length(etds)
            X1endSCS = (Y1opt(i,:)' + [mubar 0 0 0 0 0]'); %Aphelion in normalised Sun-Centered Synodic frame 
            X1endSCS = [X1endSCS(1:3)*LU ; X1endSCS(4:6)*LU/TU]; %Aphelion in SCS frame dimensionalised
            STM = cspice_sxform('SUNEARTH','ECLIPJ2000', etds(i)); %State Transformation Matrix
            XLeg1(:,i) = STM*X1endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        end
        XLeg1 = XLeg1';
        %Rotating Y2 trajectory into ECLJ2K frame
        SOIDepetd = cspice_str2et(cellstr(SOIDepDate));
        etds = linspace(EPFDepetd,SOIDepetd,length(Y2opt));
        %Rotating Y1 trajectory into ECLJ2K frame
        XLeg2 = zeros(6,length(Y2opt));
        for i = 1:1:length(etds)
            X1endSCS = (Y2opt(i,:)' + [mubar 0 0 0 0 0]'); %Aphelion in normalised Sun-Centered Synodic frame 
            X1endSCS = [X1endSCS(1:3)*LU ; X1endSCS(4:6)*LU/TU]; %Aphelion in SCS frame dimensionalised
            STM = cspice_sxform('SUNEARTH','ECLIPJ2000', etds(i)); %State Transformation Matrix
            XLeg2(:,i) = STM*X1endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        end
        XLeg2 = XLeg2';
        

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
        plot3(XLeg2(:,1), XLeg2(:,2), XLeg2(:,3),'-k','LineWidth', 3.5);
        plot3(XLeg2(end,1), XLeg2(end,2), XLeg2(end,3), 'ok', 'MarkerFaceColor', 'y','MarkerSize',8);
        %Plot transfer for Leg 3
        plot3(XLeg3(:,1), XLeg3(:,2), XLeg3(:,3),'-k','LineWidth', 3.5);
        plot3(XLeg3(end,1), XLeg3(end,2), XLeg3(end,3), 'ok', 'MarkerFaceColor', 'r','MarkerSize',8);
        %Plot orbit for asteroid 
        plot3(RV_Ast(1,:), RV_Ast(2,:), RV_Ast(3,:),'-','Color', 'b','LineWidth', 3);
        
        axis equal 
        ax = gca;
        ax.FontSize = 14; 
        ax.FontSize = 14; 
        
        title('Mission Scenario 1','FontSize',18)
        legend('','','','SCI Frame','SEL2','Manifold Transfer','Perigee','SOI Escape','SOI','NEO Transfer','Rendezvous','NEO Orbit','FontSize',16);
        exportgraphics(fig1,'D22016HO3FINALFULLLEGEND.png','Resolution',1000)
end



% % Set last_best to best
% if optimvalues.iteration == 0
% last_best = best;
% else
%     plot3(Y2opt(:,1), Y2opt(:,2), Y2opt(:,3),':','Color', 'b','LineWidth', 1.5);
% end
