function [X,DV1,DV2,DV3,TOF1,TOF2,TOF3,TOF4] = InitialGuess(Design,month,trajmin,X0min,Xfmin,Vdep,Varr,JD_Dep,X1,t1e,Xsc,RV_Ast,TOFhyp,XPerigeeMin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates Intitial Guess for Single-Shooting Method Optimisation. 
%
% INPUTS:
%   - State space at SEL2, X0 (6 x 1)
%   - Time of SEL2 escape, TOF1 (1 x 1)
%   - Departure velocity, DV1 (3 x 1)
%   - Time of transfer, TOF2 (1 x 1)
%   - Relative velocity, DV2 (3 x 1)
% OUTPUTS:                                               
% - Optimisation vector X (14 x 1), including:
%                                                                                
% Author: Sho Wright (sw01745)                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adding Paths
addpath(genpath('AsteroidSearchResults'));
addpath(genpath('Spice\mice\lib'));
addpath(genpath('Spice\mice\src\mice'));
addpath(genpath('Spice\generic_kernels'));
cspice_furnsh(which('naif0012.tls'));
cspice_furnsh(which('pck00010.tpc'));
cspice_furnsh(which('gm_de431.tpc'));
cspice_furnsh(which('de430.bsp'));
cspice_furnsh(which('SUNEARTH.txt'));

%% Defining Variables
LU = 1.4959787070e08;
muSun = 1.32712440018e11; 
n = sqrt((muSun)/LU^3); 
TU = 1/n; %Normalising Units

dep = trajmin(2); arr = trajmin(1); %Row and Column of minimum departure 
Xast = RV_Ast(:,arr); %Asteroid state space at arrival 
vast = Xast(4:6)/LU*TU; %Asteroid velocity NU
X0 = X1(1,:)'; %SEL2 state vector
TOF1 = t1e(end); %TOF Leg 1, manifold transfer
Xsc = Xsc(dep,:)'; %S/C state space at departure
vsc = Xsc(4:6);
DV1 = (X0min(4:6) - vsc)/LU*TU;
DepJD = JD_Dep(dep);
switch Design
    case 1
        TOF2 = trajmin(6)/TU; %Time to transfer to asteroid
        DV2 = (vast - Xfmin(4:6)); %Arrival V-inf
        X = [X0 ; TOF1 ; DV1 ; TOF2 ; DV2 ; DepJD]; %Total optimisation vector
        DV3 = NaN; TOF3 = NaN; TOF4 = NaN;
    case 2
        DV1 = Vdep/LU*TU;
        TOF2 = TOFhyp/TU; %Time taken to reach Earth SOI
        DV2 = zeros(3,1); %Micro-adjustments for optimiser
        TOF3 = trajmin(6)/TU - TOF2; %Time to transfer from Earth-NEO
        DV3 = Varr/LU*TU;
%         DV3 = (vast - Xfmin(4:6)); %Arrival V-inf
        X = [X0 ; TOF1 ; DV1 ; TOF2 ; DV2 ; TOF3 ; DV3 ; DepJD]; %Total optimisation vector
        TOF4 = NaN;
    case 4
        TOF2 = XPerigeeMin(month).TOF; %Time taken from LGA to EPF
        TOF3 = TOFhyp/TU; %Time taken to reach Earth SOI
        DV2 = zeros(3,1); %Micro-adjustments for optimiser
        TOF4 = trajmin(6)/TU - TOF2; %Time to transfer from Earth-NEO
        DV3 = (vast - Xfmin(4:6)); %Arrival V-inf
        X = [X0 ; TOF1 ; TOF2 ; DV1 ; TOF3 ; DV2 ; TOF4 ; DV3 ; DepJD]; %Total optimisation vector
end