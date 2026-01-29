function RV = COE2RV(coe,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 2). 
% ii). Convert the state of the satellite from COE to SCI position and 
% velocity components (R & V). Also calls upon 'Kepler' to convert Mean
% Anomaly to True anomaly.
%                                                         
% INPUTS:                                                
% - Gravitational Parameter, mu
%   - Sun mu = 1.32712440018e11 (km^3/s^2)
% - 2 Orbital Element sets:
%   - [a1,a2] = Semi-major Axis (km)
%   - [e1,e2] = Eccentricity (-)
%   - [i1,i2] = Inclination (deg)
%   - [Node1, Node2] = Longitude of Ascending Node (deg)
%   - [ArgP1,ArgP2] = Argument of Perihelion (deg)
%   - [MA1,MA2] = Mean anomaly (deg)
%                                                       
% OUTPUTS:                                               
% - Position Components: [x y z] (km)
% - Velocity Components: [Vx Vy Vz] (Km/s) 
%                                                                                 
% Author: Sho Wright (sw01745)                           
% Date created: 05/10/22                                         
% Date updated: 07/07/23                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining variables
a = coe(1,:); e = coe(2,:); i= coe(3,:); 
Node = coe(4,:); ArgP = coe(5,:); MA = coe(6,:); %Indexing Classical Orbital Elements

TA = zeros(length(MA),1);
for k = 1:length(MA)
    TA(k,:) = Kepler(e(k,:),deg2rad(MA(k,:)));
end

%% Calculating Parameters
p = a*(1-e ^2); %Semi-latus rectum
r = p/(1 + e*cosd(TA)); %Orbit equation (km)
h = sqrt(mu*a*(1-e^2));  %Specific Angular Momentum (km^2/s)
%% Coordinates in the perifocal reference system Oxyz 
% Position vector coordinates
x = r*cosd(TA);
y = r*sind(TA);
z = 0;
Perifocal_pos = [x;y;z];
% Velocity vector coordinates
Vx = -(mu/h)*sind(TA);
Vy = (mu/h)*(cosd(TA) + e);
Vz = 0;     
Perifocal_vel = [Vx;Vy;Vz];
%% Directional Cosine Matrices (DCMs)
R3_Node = [cosd(Node) sind(Node) 0 ; -sind(Node) cosd(Node) 0 ; 0 0 1]; %Rotation about 3rd axis
R1_i= [1 0 0 ; 0 cosd(i) sind(i) ; 0 -sind(i) cosd(i)]; %Rotation about 2nd axis
R3_ArgP = [cosd(ArgP) sind(ArgP) 0 ; -sind(ArgP) cosd(ArgP) 0 ; 0 0 1]; %Rotation about 1st axis
%% 3-1-3 Rotational Sequence
PJ = R3_ArgP*R1_i*R3_Node;
JP = transpose(PJ);
%% Position and velocity vector in SCI
SCI_pos = JP*Perifocal_pos;
SCI_vel = JP*Perifocal_vel;
RV = [SCI_pos ; SCI_vel]; %Defining X to be called upon by 'OrbitalPropagation.m'
end