function [c,ceq] = nonlcon(X,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Completes the single shooting method for the use of fmincon for leg 1 
% and leg 2 of design 1.
%
% INPUTS:
% - Optimisation vector X (14 x 1), including:
%   - State space at SEL2, X0 (6 x 1)
%   - Time of SEL2 escape, TOF1 (1 x 1)
%   - Departure velocity, DV1 (3 x 1)
%   - Time of transfer, TOF2 (1 x 1)
%   - Relative velocity, DV2 (3 x 1)
% OUTPUTS:                                               
% - Nonlinear equality constraints, ceq
% - Nonlinear inequality constraints, c
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 08/08/23                                         
% Date updated: 15/10/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unpacking Parameters
Design = param.Design;
sidereal = param.sidereal;          
LU = param.LU; TU = param.TU;
mubar = param.mubar; mu = param.mu*TU^2/LU^3;    
RV_Ast = param.RV_Ast;scenario = param.scenario; 
JD_Arr = param.JD_Arr; TOF1ori = param.TOF1*TU;
Event = param.Event; TOF2ori = param.TOF2*TU; 
switch Design
    case 2
        TOF3ori = param.TOF3;
    case 4
        TOF3ori = param.TOF3;
end

%% CR3BP Domain Transfers
%Leg 1 - Escaping SEL2 (All designs)
options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16); %Setting tolerances (varied for situation)
[TOF1opt,Y1opt] = ode45(@(t,X)CR3BPEOM(mubar,X,t),[0 X(7)],X(1:6),options); %Integrating system of 1st order ODEs
TOF1opt = TOF1opt(end)*TU; %Dimensionalise
ropt1 = Y1opt(end,1:3)'; vopt1 = Y1opt(end,4:6)'; %For Aphelion/Perigee constraints

switch Design
    case 2 %Rotating DV1 into Synodic, CR3BP EOM to SOI Edge
        %Setting DV1 in Ecliptic J2K to a "state"
        DV1EJ2K = [zeros(3,1)*LU; X(8:10)*LU/TU];
        %Using Spice syntax
%         Time = X(19)*sidereal - TOF1ori + TOF1opt - 2451545.0*sidereal;
        etd1 = X(19)*sidereal - TOF1ori + TOF1opt;% - 2451545.0*sidereal;
        DateDep = datetime(etd1/sidereal,'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');
        etd1 = cspice_str2et(cellstr(DateDep));
        STM = cspice_sxform('ECLIPJ2000','SUNEARTH', etd1); %State Transformation Matrix 
%         RV_Earth = cspice_spkezr('399',Time,'ECLIPJ2000','NONE','10');
%         [r,v] = Rotation(DV1EJ2K',RV_Earth',mubar,'CR3BP-SCI');
%         DV1Syn = [r;v];
        DV1Syn = DV1EJ2K'*STM; %DV1 in Synodic frame
        DV1Syn = [DV1Syn(1:3)/LU , DV1Syn(4:6)/LU*TU];
        Y1new = Y1opt(end,:)' + DV1Syn';
        [TOF2opt,Y2opt] = ode45(@(t,X)CR3BPEOM(mubar,X,t),[0 X(11)],Y1new,options); %Integrating system of 1st order ODEs
        ropt2 = Y2opt(end,1:3)'; vopt2 = Y2opt(end,4:6)'; %For SOI constraints
        TOF2opt = TOF2opt(end)*TU; 
    case 4 
        %Leg 2 - LGA at Perilune to Perigee
        [TOF2opt,Y2opt] = ode45(@(t,X)CR3BPEOM(mubar,X,t),[0 X(8)],Y1opt(end,:)',options); %Integrating system of 1st order ODEs
        TOF2opt = TOF2opt(end)*TU;
        ropt2 = Y2opt(end,1:3)'; vopt2 = Y2opt(end,4:6)'; %For Perigee constraints
        %Setting DV1 in Ecliptic J2K to a "state"
        DV1EJ2K = [zeros(3,1) ; X(9:11)];
        %Using Spice syntax
        Time = X(20)*sidereal - TOF1ori - TOF2ori + TOF1opt + TOF2opt - 2451545.0*sidereal;
        STM = cspice_sxform('ECLIPJ2000','SUNEARTH', Time); %State Transformation Matrix  
        DV1Syn = DV1EJ2K*STM; %DV1 in Synodic frame
        Y2new = Y2opt(end,:)' + DV1Syn;
        [TOF3opt,Y3opt] = ode45(@(t,X)CR3BPEOM(mubar,X,t),[0 X(12)],Y2new,options); %Integrating system of 1st order ODEs
        ropt3 = Y3opt(end,1:3)'; %For SOI constraints
        TOF3opt = TOF3opt(end)*TU; 
end

%% CR3BP - 2BP Domain Rotation
switch Design 
    case 1 
        %Using Spice syntax
        Time = X(15)*sidereal - TOF1ori + TOF1opt;
        Time = Time - 2451545.0*sidereal; %2451545.0 = JD 1st Jan 
        %Rotating end of Leg 1 (Aphelion) into SCI frame
        Y1endSCS = (Y1opt(end,:)' + [mubar 0 0 0 0 0]'); %Aphelion in normalised Sun-Centered Synodic frame 
        Y1endSCS = [Y1endSCS(1:3)*LU ; Y1endSCS(4:6)*LU/TU];
        STM = cspice_sxform('SUNEARTH','ECLIPJ2000', Time); %State Transformation Matrix  
        Y1 = STM*Y1endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        Y1 = [Y1(1:3)/LU ; Y1(4:6)/LU*TU];
    case 2
        %Using Spice syntax
        Time = X(19)*sidereal - TOF1ori + TOF1opt - TOF2ori + TOF2opt;% - 2451545.0*sidereal;
        Time = Time - 2451545.0*sidereal;
%         DateDep = datetime(etd2/sidereal,'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');
%         etd2 = cspice_str2et(cellstr(DateDep));
        %Rotating end of Leg 3 (SOI) into SCI frame
        Y2endSCS = (Y2opt(end,:)' + [mubar 0 0 0 0 0]'); %SOI Edge in normalised Sun-Centered Synodic frame 
        Y2endSCS = [Y2endSCS(1:3)*LU ; Y2endSCS(4:6)*LU/TU];
        STM = cspice_sxform('SUNEARTH','ECLIPJ2000', Time); %State Transformation Matrix  
        Y2 = STM*Y2endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        Y2 = [Y2(1:3)/LU ; Y2(4:6)/LU*TU];
    case 4
        %Using Spice syntax
        Time = Time - TOF3ori + TOF3opt - 2451545.0*sidereal;
        %Rotating end of Leg 4 (SOI) into SCI frame
        Y3endSCS = (Y3opt(end,:)' + [mubar 0 0 0 0 0]'); %SOI Edge in normalised Sun-Centered Synodic frame 
        Y3endSCS = [Y3endSCS(1:3)*LU ; Y3endSCS(4:6)*LU/TU];
        STM = cspice_sxform('SUNEARTH','ECLIPJ2000', Time); %State Transformation Matrix  
        Y3 = STM*Y3endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        Y3 = [Y3(1:3)/LU ; Y3(4:6)/LU*TU];
end

%% 2BP Domain Transfers
switch Design 
    case 1
        X1 = Y1 + [zeros(3,1) ; X(8:10)]; %Transfer state space vector
        [~, Y2opt] = ode45(@(t,X)EOM2BP(t,X,mu),[0 X(11)],X1,options); %Integrating system of 1st order ODEs
        Yend = Y2opt(end,:)';
        JDcon = (Time + X(11)*TU)/sidereal + 2451545.0;
    case 2
        X2 = Y2 + [zeros(3,1) ; X(12:14)]; %Transfer state space vector
        [~, Y3opt] = ode45(@(t,X)EOM2BP(t,X,mu),[0 X(15)],X2,options); %Integrating system of 1st order ODEs
        Yend = Y3opt(end,:)';
        JDcon = (Time + X(15)*TU)/sidereal + 2451545.0;
    case 4
        X3 = Y3 + [zeros(3,1) ; X(13:15)]; %Transfer state space vector
        [~, Y4opt] = ode45(@(t,X)EOM2BP(t,X,mu),[0 X(16)],X3,options); %Integrating system of 1st order ODEs
        Yend = Y4opt(end,:)';
        JDcon = (Time + X(16)*TU)/sidereal + 2451545.0;
end
RV_Ast = spline(JD_Arr,RV_Ast,JDcon);
RV_Ast = [RV_Ast(1:3)/LU ; RV_Ast(4:6)*TU/LU]; %Normalised Asteroid State Space 

%% Defining Nonlinear Constraints
switch Design
    case 1
         Y1new = NaN; Y2new = NaN; ropt2 = NaN; vopt2 = NaN; ropt3 = NaN;
        [c,ceq] = constraintsinput(Design,X,scenario,Event,Yend,Y1new,Y2new,RV_Ast,ropt1,vopt1,ropt2,vopt2,ropt3,mubar,TU,LU,mu);
    case 2
        Y2new = NaN; ropt3 = NaN;
        [c,ceq] = constraintsinput(Design,X,scenario,Event,Yend,Y1new,Y2new,RV_Ast,ropt1,vopt1,ropt2,vopt2,ropt3,mubar,TU,LU,mu);
    case 4
        [c,ceq] = constraintsinput(Design,X,scenario,Event,Yend,Y1new,Y2new,RV_Ast,ropt1,vopt1,ropt2,vopt2,ropt3,mubar,TU,LU,mu);
end