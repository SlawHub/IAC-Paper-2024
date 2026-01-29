function MOIDMain(mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                          
% Returns a list of MOID of each orbit in a given asteroid database using 
% orbital elements. 
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
% - List of asteroids with suitable MOID (au)                            
%                                                        
% Author: Sho Wright (sw01745)                           
% Date created: 07/07/23                                         
% Date updated: 07/07/23                                                                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting up script
close all; clc; format short;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 16);
%% Defining Constants
AU = 1.4959787070e08;
n = sqrt(mu/AU^3); mubar = 3.0032080443e-06;
c = 299792; % in km/s
%% Retrieving Earth State
fid = fopen('Earth.txt');
A = textscan(fid,'%f %s %f %f %f %f %f %f', 'Delimiter', ',')';
fclose(fid);
% JD_Dep = A{1}';
% Datedep = A{2}';
RV_Earth = [A{3}, A{4}, A{5} A{6}, A{7}, A{8}];

choice1 = input('Enter 1 S/C design 1, 2 for S/C design 2, or 3 for Earth COE: '); %User input configuration selection

choice2 = input('Enter 1 for JPL Data, 2 for MPC Data: '); %User input configuration selection

%% Retrieving Spacecraft Data
switch choice1
    case 1
        load("CR3BPOutput1.mat",'Ye')
    case 2
        load("CR3BPOutput2.mat",'Ye')
end

%Calculating MOID threshold
threshold = MOIDthreshold(RV_Earth,Ye,mubar);
threshold = round(threshold/AU,3);

switch choice1 
    case {1,2}   
        % Rotating Spacecraft into Inertial Frame
        r_sci = zeros(length(RV_Earth),3);
        v_sci = zeros(length(RV_Earth),3);
        SCcoe = zeros(length(RV_Earth),6);
    
        for i = 1:length(RV_Earth)
            [r_sci(i,:),v_sci(i,:)] = Rotation(Ye(end,:),RV_Earth(i,:),mubar); %Rotates starting position at Perihelion into inertial frame
            r_sci(i,:) = r_sci(i,:)*AU; v_sci(i,:) = v_sci(i,:)*(AU*n);
            SCcoe(i,:)= RV2COE(r_sci(i,:)', v_sci(i,:)',mu);
        end

        e1 = SCcoe(:,2); V1 = SCcoe(:,6);
        MA1 = atan2d(-sqrt(1-e1.^2).*sind(V1),-e1-cosd(V1)) + 180 - e1.*sqrt(1-e1.^2).*sind(V1)./1+cosd(V1);
        SCcoe = [SCcoe(:,1:5)' ; MA1'];
    case 3
        %% Retrieving Earth Data
        fid = fopen('EarthCOE.txt');
        A = textscan(fid,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f', 'Delimiter', ',')';
        fclose(fid);
        a1 = A{12}; e1 = A{3}; i1 = A{5}; Node1 = A{6}; ArgP1 = A{7}; MA1 = A{10};
        Earthcoe = [a1(1); e1(1); i1(1); Node1(1); ArgP1(1); MA1(1)];
end

switch choice2
    case 1
        %% Retrieving Asteroid Data - JPL Horizons
        A = readtable('sbdb_query_results_27_07_23.csv',"VariableNamingRule","preserve");
        anomaly = find(strcmp('       (2002 PD153)',A{:,2}));
        A(anomaly,:) = []; %#ok<FNDSB> 
        a2 = table2array(A(:,34))*AU; 
        e2 = table2array(A(:,33)); 
        i2 = table2array(A(:,36));  
        Node2 = table2array(A(:,37)); 
        ArgP2 = table2array(A(:,38)); 
        MA2 = table2array(A(:,39));
        %EarthMOID = table2array(A(:,10));
        AstName = table2array(A(:,2));
        %NoSat = table2array(A(:,8));
        AbsMag= table2array(A(:,9));
        albedo = table2array(A(:,18));
        %Diam = table2array(A(:,16));
        Astcoe = [a2'; e2'; i2'; Node2'; ArgP2'; MA2'];
    case 2
        %% Retrieving Asteroid Data - MPC
        %load("A.mat","A");
        A = readtable('MPCORB_27_07_23.DAT','FileType','fixedwidth');
        a2 = table2array(A(:,11))*AU; 
        e2 = table2array(A(:,9)); 
        i2 = table2array(A(:,8));  
        Node2 = table2array(A(:,7)); 
        ArgP2 = table2array(A(:,6)); 
        MA2 = table2array(A(:,5));
        AstName = table2array(A(:,23));
        AbsMag = table2array(A(:,2));
        %Pulling out NEOs
        rp = a2/AU.*(1-e2); ra = a2/AU.*(1+e2); index = 1:1:length(a2);
        NEO = [rp, ra, index'];
        %NEO = NEO((NEO(:,2)>0.85),:);
        NEO = NEO((NEO(:,1)<1.3),:);
        
        %Resizing arrays
        a2 = a2(NEO(:,3),:); e2 = e2(NEO(:,3),:); i2 = i2(NEO(:,3),:);
        Node2 = Node2(NEO(:,3),:); ArgP2 = ArgP2(NEO(:,3),:); MA2 = MA2(NEO(:,3));
        AstName = AstName(NEO(:,3),:); AbsMag = AbsMag(NEO(:,3),:);
        Astcoe = [a2'; e2'; i2'; Node2'; ArgP2'; MA2'];
end
% Anom = A(find(e2(:,:) == 0),:);
%% Running MOID Calculator
%Initialising
RV1 = zeros(6,length(Astcoe));
RV2 = zeros(6,length(Astcoe));
MOID = zeros(length(Astcoe),1);
Vsc = zeros(length(Astcoe),1);
Vast = zeros(length(Astcoe),1);
Vdot = zeros(length(Astcoe),1);
Vrel = zeros(length(Astcoe),1);

for i = 1:1:length(Astcoe)
    if i2(i,:) < 10 && e2(i,:) < 1 
        switch choice1 
            case {1,2}
                [MOID(i,:),RV1(:,i),RV2(:,i)] = calculate_moid(SCcoe(:,1),Astcoe(:,i),mu);
                Vsc(i,:) = norm(RV1(4:6,i))'; Vast(i,:) = norm(RV2(4:6,i))'; Vdot(i,:) = dot(RV1(4:6,i)',RV2(4:6,i)');
                Vrel(i,:) = dot(sqrt(1 - (c^2 - Vsc(i,:)^2)*(c^2 - Vast(i,:)^2)/(c^2 -Vdot(i,:))^2),c);
            case 3
                [MOID(i,:),RV1(:,i),RV2(:,i)] = calculate_moid(Earthcoe(:,1),Astcoe(:,i),mu);
        end
    else
        MOID(i,:) = NaN;
        Vsc(i,:) = NaN; Vast(i,:) = NaN; 
        Vdot(i,:) = NaN; Vrel(i,:) = NaN;
        continue
    end
end

%% Outputting Final List
% Vrel = vecnorm(RV1(4:6,:) - RV2(4:6,:));
index = 1:1:length(Astcoe);
MOID = [MOID/AU, index'];
MOID = MOID((MOID(:,1)<threshold),:);
MOID = MOID((MOID(:,1)>0),:);
Astcoe = [Astcoe(1,:)/AU ; Astcoe(2:6,:)];

switch choice2
    case 1
        AstDat = [num2cell(MOID(:,2)), AstName(MOID(:,2),:),num2cell(MOID(:,1)),num2cell(AbsMag(MOID(:,2),:)),num2cell(albedo(MOID(:,2),:)),num2cell(Vrel(MOID(:,2),:)), num2cell(Astcoe(:,MOID(:,2))')];
        tnames = ["Index","Asteroid" ,"MOID (AU)","Absolute Magnitude","Albedo","Relative velocity (km/s)", "a (AU)", "e (-)", "i (deg)", "LOAN (deg)", "ArgP (deg)","MA (deg)"]; %"RV1","RV2    
    case 2
        AstDat = [num2cell(MOID(:,2)), AstName(MOID(:,2),:), num2cell(MOID(:,1)),num2cell(AbsMag(MOID(:,2),:)),num2cell(Vrel(MOID(:,2),:)),num2cell(Astcoe(:,MOID(:,2))')];
        tnames = ["Index","Asteroid" ,"MOID (AU)", "Absolute Magnitude","Relative velocity (km/s)", "a (AU)", "e (-)", "i (deg)", "LOAN (deg)", "ArgP (deg)","MA (deg)"]; %"RV1","RV2
end

AstTable = array2table(AstDat,"VariableNames",tnames);
disp(AstTable)
fprintf('Returning %d asteroids \n',height(AstTable))