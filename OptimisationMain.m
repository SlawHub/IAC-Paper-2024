function Xopt = OptimisationMain(spkid,AstName,Design,month)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements constrained direct single shooting method using fmincon to 
% minimise the objective functions delta-v and TOF for leg 1 and 2 of D1.
%
% INPUTS:
% - Optimisation vector X (14 x 1), including:
%   - State space at SEL2, X0 (6 x 1)
%   - Time of SEL2 escape, TOF1 (1 x 1)
%   - Departure velocity, DV1 (3 x 1)
%   - Time of transfer, TOF2 (1 x 1)
%   - Relative velocity, DV2 (3 x 1)
% OUTPUTS:                                               
% - Optimised flight times, TOF1 and TOF2 
% - Optimised departure and arrival delta-vs, DV1 and DV2 
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 08/08/23                                         
% Date updated: 15/10/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting up script
close all; clc; format short;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 16);

%% Defining Variables
mubar = 3.0032080443e-06;
sidereal = 86400; LU = 1.4959787070e08;
muSun = 1.32712440018e11; n = sqrt((muSun)/LU^3); 
JDJ2000 = 2451545.0;
TU = 1/n; %Normalising Units

%% Adding NAIF Spice Path
addpath(genpath('Spice\mice\lib'));
addpath(genpath('Spice\mice\src\mice'));
addpath(genpath('Spice\generic_kernels'));
cspice_furnsh(which('naif0012.tls'));
cspice_furnsh(which('pck00010.tpc'));
cspice_furnsh(which('gm_de431.tpc'));
cspice_furnsh(which('de430.bsp'));
cspice_furnsh(which('SUNEARTH.txt'));

%% Retrieving Spacecraft Data
switch Design
    case 1
        load('D1CR3BP.mat','X1','t1e','X1end');
    case 2
        load('D2CR3BP.mat','X1','X1end','t1e')
%         load('possible_targets.mat','sc_state')
%         X1 = sc_state';
    case 4
        load('D4CR3BP.mat','X1','X1end','t1e','XPerigeeMin')
        load('DateFinder.mat','AlignedTable')
end

%% Retrieving Earth Data
%Earth state in day resolution
fid = fopen('Earth.txt'); 
A = textscan(fid,'%f %s %f %f %f %f %f %f', 'Delimiter', ',')';
fclose(fid);
JD_Dep = A{1}';
DATEdep = A{2}';
RV_Earth = [A{3} A{4} A{5} A{6} A{7} A{8}];
%Slicing out tick labels for each axis
DATEdep = cellfun(@(x)x(1:end-13),DATEdep,'un',0);
DATEdep = cellfun(@(x)x(6:end),DATEdep,'un',0);

%% Retrieving Asteroid State
rmpath(genpath('AsteroidSearchResults'))
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

%% Retrieving Transfer Data
transfername = sprintf('D%dTransfer%s.mat',Design,AstName);
switch Design 
    case 1
        load(transfername,'X0min','Xfmin','trajmin');
    case 2
        load(transfername,'X0min','Xfmin','trajmin','Vdep','Varr'); 
end
dep = trajmin(2,:); arr = trajmin(1,:);

%% Generating S/C Orbit
switch Design
    case 1
        X1endSCS = (X1end(end,:)' + [mubar 0 0 0 0 0]'); %Aphelion in normalised Sun-Centered Synodic frame 
        X1endSCS = [X1endSCS(1:3)*LU ; X1endSCS(4:6)*LU/TU]; %Aphelion in SCS frame dimensionalised
        Xsc = zeros(6,length(DATEdep)); %Initialising
        etd = zeros(length(DATEdep),1);
        for i = 1:length(DATEdep)
            etd(i,:) = cspice_str2et(DATEdep(i)); %Ephemeris time at Aphelion
            STM = cspice_sxform('SUNEARTH','ECLIPJ2000', etd(i,:)); %State Transformation Matrix
            Xsc(:,i) = STM*X1endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        end
        Xsc = Xsc';
    case 2
        X1endSCS = (X1end(end,:)' - [1-mubar 0 0 0 0 0]'); %Aphelion in normalised Sun-Centered Synodic frame 
        X1endSCS = [X1endSCS(1:3)*LU ; X1endSCS(4:6)*LU/TU]; %Aphelion in SCS frame dimensionalised
        Xsc = zeros(6,length(DATEdep)); %Initialising
        etd = cspice_str2et(DATEdep); %Ephemeris time at Aphelion
        etd = etd';
        for i = 1:length(DATEdep)
            STM = cspice_sxform('SUNEARTH','ECLIPJ2000', etd(i,:)); %State Transformation Matrix
            Xsc(:,i) = STM*X1endSCS; %Rotating into SCI (Ecliptic J2000 frame)
        end
        Xsc = Xsc';
        etPerigee = trajmin(8);
        TOFhyp = etd(dep) - etPerigee;
    case 4
        X2end = XPerigeeMin(month).Xleg2(:,end);
        X2endSCS = (X2end(end,:)' - [1-mubar 0 0 0 0 0]'); %Aphelion in normalised Sun-Centered Synodic frame 
        X2endSCS = [X2endSCS(1:3)*LU ; X2endSCS(4:6)*LU/TU]; %Aphelion in SCS frame dimensionalised
        LGADate = datetime(string(table2cell(AlignedTable(month,3))));
        Tleg2 = XPerigeeMin(month).TOF*TU;
        EPFDate = LGADate + duration(seconds(Tleg2));
        etd = cspice_str2et(cellstr(EPFDate)); %Ephemeris time at Aphelion
        STM = cspice_sxform('SUNEARTH','ECLIPJ2000', etd); %State Transformation Matrix
        Xsc = STM*X2endSCS; %Rotating into SCI (Ecliptic J2000 frame)
end

%% Building Optimisation vector 
switch Design
    case 1
        month = NaN; TOFhyp = NaN; XPerigeeMin = NaN; Vdep = NaN; Varr = NaN;
        [X,DV1,DV2,~,TOF1,TOF2,~,~] = InitialGuess(Design,month,trajmin,X0min,Xfmin,Vdep,Varr,JD_Dep,X1,t1e,Xsc,RV_Ast,TOFhyp,XPerigeeMin);
    case 2
        month = NaN; XPerigeeMin = NaN;
        [X,DV1,DV2,DV3,TOF1,TOF2,TOF3,~] = InitialGuess(Design,month,trajmin,X0min,Xfmin,Vdep,Varr,JD_Dep,X1,t1e,Xsc,RV_Ast,TOFhyp,XPerigeeMin);      
    case 4
        [X,DV1,DV2,DV3,TOF1,TOF2,TOF3,TOF4] = InitialGuess(Design,month,trajmin,X0min,Xfmin,Vdep,Varr,JD_Dep,X1,t1e,Xsc,RV_Ast,TOFhyp,XPerigeeMin);
end

%% Defining Linear Constraints
Adiag = diag(ones(6,1));
Azeros = zeros(6,length(X)-6);
%Linear equality constraints
Aeq = horzcat(Adiag,Azeros);
beq = X(1:6); 
%Linear inequality constraints
A = []; b = [];

%% Defining Boundary Conditions
switch Design
    case 1    
        lb = -Inf*ones(15,1); %Lower Bounds
        lb(7) = 0; %Setting TOF1 = 0
        lb(11) = 0; %Setting TOF2 = 0
        lb(15) = JD_Dep(dep-1);
        ub = Inf*ones(15,1); %Upper Bounds
        ub(15) = JD_Dep(dep+1);
    case 2
        lb = -Inf*ones(19,1); %Lower Bounds
        lb(7) = 0; %Setting TOF1 = 0
        lb(11) = 0; %Setting TOF2 = 0
        lb(15) = 0; %Setting TOF3 = 0
        lb(19) = JD_Dep(dep-1); %JD minus 1 day
        ub = Inf*ones(19,1); %Upper Bounds
        ub(19) = JD_Dep(dep+1); %JD plus 1 day
    case 4
        lb = -Inf*ones(20,1); %Lower Bounds
        lb(7) = 0; %Setting TOF1 = 0
        lb(8) = 0;%Setting TOF2 = 0
        lb(12) = 0; %Setting TOF3 = 0
        lb(16) = 0; %Setting TOF4 = 0
        lb(20) = JD_Dep(dep-1); %JD minus 1 day
        ub = Inf*ones(19,1); %Upper Bounds
        ub(20) = JD_Dep(dep+1); %JD plus 1 day
end

%% Defining original objectives
DateDepi = datetime(trajmin(5),'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');
DateArri = datetime(trajmin(4),'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');
TOF1i = TOF1*TU/sidereal;
DV1i = norm(DV1)*LU/TU;
TOF2i = TOF2*TU/sidereal;
DV2i = norm(DV2)*LU/TU;
switch Design
    case 1
        totTOFi = TOF1i + TOF2i;
        totDVi = DV1i;
    case 2
        DV3i = norm(DV3)*LU/TU;
        TOF3i = TOF3*TU/sidereal;
        totTOFi = TOF1i + TOF2i + TOF3i;
        totDVi = DV1i + DV2i;
    case 4
        DV3i = norm(DV3)*LU/TU;
        TOF3i = TOF3*TU/sidereal;
        TOF4i = TOF4*TU/sidereal;
        totTOFi = TOF1i + TOF2i + TOF3i+TOF4i;
        totDVi = DV1i + DV2i;
end

%% Calling fmincon
optimise = 1;
I = 1;
while optimise == 1
    Event = input(['Keep Orbit Event Function Conditions?: \n' ...
                   '1 for Yes, 2 for No: ']);

    %Cost Function selection
    cost(I) = input(['Choose Cost Function: \n' ...
                     '1. Total Departure Velocity: \n' ...
                     '2. Total TOF: \n' ...
                     '3. Total Delta-V: \n' ...
                     '4. Arrival Velocity: \n']); %#ok<AGROW>
    switch Design %Cost Function selection storage 
        case 1
            if cost(I) == 1
                costfun = @(X) norm(X(8:10)); %1e-10^2* regularisation term
                coststr(I) = cellstr("DV1"); %#ok<AGROW> 
            elseif cost(I) == 2
                costfun = @(X) (X(7) + X(11)); 
                coststr(I) = cellstr("TOTAL TOF"); %#ok<AGROW> 
            elseif cost(I) == 3
                costfun = @(X) norm(X(8:10))+ norm(X(12:14));
                coststr(I) = cellstr("TOTAL DV"); %#ok<AGROW> 
            elseif cost(I) == 4
                costfun = @(X) norm(X(12:14));
                coststr(I) = cellstr("DV2"); %#ok<AGROW> 
            end
        case 2
           if cost(I) == 1
                costfun = @(X) norm(X(8:10))+ norm(X(12:14)); %1e-10^2* regularisation term
                coststr(I) = cellstr("DV1"); %#ok<AGROW> 
            elseif cost(I) == 2
                costfun = @(X) (X(7) + X(11) + X(15));
                coststr(I) = cellstr("TOTAL TOF"); %#ok<AGROW> 
            elseif cost(I) == 3
                costfun = @(X) norm(X(8:10))+ norm(X(12:14)) + norm(X(16:18));
                coststr(I) = cellstr("TOTAL DV"); %#ok<AGROW> 
           elseif cost(I) == 4
                costfun = @(X) norm(X(16:18));
                coststr(I) = cellstr("DV3"); %#ok<AGROW>  
           end
        case 4
           if cost(I) == 1
                costfun = @(X) norm(X(9:11))+ norm(X(13:15)); %1e-10^2* regularisation term
                coststr(I) = cellstr("DV1"); %#ok<AGROW> 
            elseif cost(I) == 2
                costfun = @(X) (X(7) + X(8) + X(12) + X(16));
                coststr(I) = cellstr("TOTAL TOF"); %#ok<AGROW> 
            elseif cost(I) == 3
                costfun = @(X) norm(X(9:11))+ norm(X(13:15)) + norm(X(17:19));
                coststr(I) = cellstr("TOTAL DV"); %#ok<AGROW> 
           elseif cost(I) == 4
                costfun = @(X) norm(X(17:19));
                coststr(I) = cellstr("DV3"); %#ok<AGROW>  
           end
    end  
    switch Design %Choose Nonlinear Constraints
        case 1
            scenario(I) = input(['Choose Nonlinear constraints to enforce: \n' ...
                                 '(S/C and Asteroid positions constraint by default) \n'...
                                 'Nonlinear Inequality Constraints (c = [] by default): \n' ...
                                 '1. Departure Velocities <= 0.15 km/s \n' ...
                                 '2. Relative Velocities <= 5km/s \n']); %#ok<AGROW> 
        case {2,4}
            scenario(I) = input(['Choose Nonlinear constraints to enforce: \n' ...
                                 '(S/C, Asteroid positions, and SOI exit constraint by default) \n'...
                                 'Nonlinear Inequality Constraints: \n' ...
                                 '1. S/C is unbounded after DV1 at Perigee \n'...
                                 '2. Departure Velocities <= 0.15 km/s \n' ...
                                 '3. Relative Velocities <= 5km/s \n']); %#ok<AGROW>  
    end

    options = optimoptions('fmincon','Display','Iter',"EnableFeasibilityMode",true,'MaxIterations',500, ...
                           'MaxFunctionEvaluations',Inf,'TolFun',1e-06,'ConstraintTolerance',1e-06,"SubproblemAlgorithm","cg");  
    Algo = input(['Choose Algorithm: \n' ...
                  '1. Interior-point \n' ...
                  '2. SQP \n']);
    switch Algo
        case 1
            options.Algorithm = 'interior-point';
        case 2
            options.Algorithm = 'sqp';
    end  
    Plot = input(['Orbit Plot?: \n' ...
                  '1 for Yes, 2 for No: ']);
    switch Plot
        case 1
            switch Design
                case 1
                    options.OutputFcn = @OrbitoutfcnD1;
                case 2
                    options.OutputFcn = @OrbitoutfcnD2;
                case 4
        %             options.OutputFcn = @OrbitoutfcnD1;
            end
    end
    %Packing up parameters
    param.Design = Design; param.sidereal = sidereal;       
    param.n = n; param.mubar = mubar;
    param.mu = muSun; param.RV_Ast = RV_Ast;  
    param.dep = dep; param.arr = arr;
    param.Xsc = Xsc; param.RV_Earth = RV_Earth;   
    param.JD_Dep = JD_Dep; param.JD_Arr = JD_Arr;
    param.TOF1 = TOF1; param.TOF2 = TOF2; param.TU = TU;
    param.LU = LU; param.scenario = scenario(I);
    param.Event = Event;
    switch Design
        case 2
            param.DV3 = DV3;
            param.TOF3 = TOF3;
        case 4
            param.DV3 = DV3;
            param.TOF3 = TOF3;
            param.TOF4 = TOF4;
    end
   
    %Calling fmincon
    [Xopt,~,flag(I),output] = fmincon(@(X)costfun(X), X, A, b, Aeq, beq, lb, ub, @(X)nonlcon(X,param),options); %#ok<AGROW> 
    %Packing up data for results table
    Rithm(I) = cellstr(output.algorithm); %#ok<AGROW> 
    iter(I) = output.iterations; %#ok<AGROW> 
    
    %Defining optimised objectives
    TOF1o(I) = (Xopt(7)*TU)/sidereal; %#ok<AGROW> 
    switch Design
        case 1
            DV1o(I) = norm(Xopt(8:10))*LU/TU; %#ok<AGROW>
            TOF2o(I) = (Xopt(11)*TU)/sidereal;%#ok<AGROW> 
            DV2o(I) = norm(Xopt(12:14))*LU/TU;%#ok<AGROW> 
            totTOFo(I) = TOF1o(I) + TOF2o(I);%#ok<AGROW>]
            totDVo(I) = DV1o(I); %#ok<AGROW> 
            DateDepo(I) = datetime(Xopt(15),'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');%#ok<AGROW>
            DateArro(I) = DateDepo(I) + duration(seconds(Xopt(11)*TU));%#ok<AGROW>
        case 2
            DV1o(I) = norm(Xopt(8:10))*LU/TU; %#ok<AGROW> 
            TOF2o(I) = (Xopt(11)*TU)/sidereal;%#ok<AGROW> 
            DV2o(I) = norm(Xopt(12:14))*LU/TU;%#ok<AGROW> 
            TOF3o(I) = (Xopt(15)*TU)/sidereal;%#ok<AGROW> 
            DV3o(I) = norm(Xopt(16:18))*LU/TU;%#ok<AGROW> 
            totTOFo(I) = TOF1o(I) + TOF2o(I) + TOF3o(I);%#ok<AGROW>
            totDVo(I) = DV1o(I) + DV2o(I); %#ok<AGROW> 
            DateDepo(I) = datetime(Xopt(19),'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');%#ok<AGROW>
            DateArro(I) = DateDepo(I) + duration(seconds((Xopt(11)+Xopt(15))*TU));%#ok<AGROW>
        case 4
            TOF2o(I) = (Xopt(8)*TU)/sidereal;%#ok<AGROW> 
            DV1o(I) = norm(Xopt(9:11))*LU/TU; %#ok<AGROW> 
            TOF3o(I) = (Xopt(12)*TU)/sidereal;%#ok<AGROW> 
            DV2o(I) = norm(Xopt(13:15))*LU/TU;%#ok<AGROW> 
            TOF4o(I) = (Xopt(16)*TU)/sidereal;%#ok<AGROW> 
            DV3o(I) = norm(Xopt(17:19))*LU/TU;%#ok<AGROW> 
            totTOFo(I) = TOF1o(I) + TOF2o(I) + TOF3o(I) + TOF4o(I);%#ok<AGROW>
            totDVo(I) = DV1o(I) + DV2o(I); %#ok<AGROW> 
            DateDepo(I) = datetime(Xopt(20),'ConvertFrom','juliandate','Format','dd-MMM-uuuu HH:mm:ss');%#ok<AGROW>
            DateArro(I) = DateDepo(I) + duration(seconds((Xopt(12)+Xopt(16))*TU));%#ok<AGROW>
    end
    %Printing this iterations results
    fprintf('Original times/delta-vs: \n')
    fprintf('Total TOF: %.2f days (~ %.2f months) \n',totTOFi,totTOFi/30)
    fprintf('Total Departure delta-v: %.4f km/s\n',totDVi)
    switch Design
        case 1
            fprintf('Total Approach delta-v: %.4f km/s\n',DV2i)
        case {2,4}
            fprintf('Total Approach delta-v: %.4f km/s\n',DV3i)
    end
    fprintf('Departure Date: %s \n',DateDepi)
    fprintf('Arrival Date: %s \n',DateArri)
    fprintf('Optimal times/delta-vs: \n')
    fprintf('Total TOF: %.2f days (~ %.2f months) \n',totTOFo(I),totTOFo(I)/30)
    fprintf('Departure delta-v: %.4f km/s\n',totDVo(I))
    switch Design
        case 1
            fprintf('Approach delta-v: %.4f km/s\n',DV2o(I))
        case {2,4}
            fprintf('Approach delta-v: %.4f km/s\n',DV3o(I))
    end
    fprintf('Departure Date: %s \n',DateDepo(I))
    fprintf('Arrival Date: %s \n',DateArro(I))

    optimise = input('Choose 1 to continue optimising, or 0 to end: \n');
    switch optimise
        case 1
            I = I + 1;
            X = Xopt;
    end
end

if I == 1
    Xopt = X;
end

%% Export Option
Export = input('Export Optimisation Table? 1 for Yes, 2 for No: ');

%% Compiling Results
empty = cell(1,1);
OptCount = num2cell(1:1:length(TOF1o));
AstName = repmat(cellstr(AstName),length(OptCount),1);
spkid = repmat(num2cell(spkid),length(OptCount),1);
switch Design
    case 1
        %Exporting Original Results
        OriDat = [cellstr(AstName(1,:)),num2cell(spkid(1,:)),empty,empty,empty,empty,empty,empty,cellstr(char(DateDepi)),cellstr(char(DateArri)),...
                  num2cell(TOF1i(:)),num2cell(TOF2i(:)),num2cell(totTOFi(:)),num2cell(totTOFi(:)/30),num2cell(DV1i(:)), num2cell(DV2i(:))];
        tnames1 = ["Name","Spkid","Opt.Iter","Cost-Fcn","Constraints","Algorithm","No.Iterations","Exit-Flag","Dep Date","Arr Date",...
                  "TOF1 (days)","TOF2 (days)","TOF total (days)","TOF total (months)","Dep.Vel (km/s)","Rel.Vel (km/s)"];
        %Exporting Optimsation Data
        OptDat = [AstName(:),spkid(:),OptCount(:),coststr(:),num2cell(scenario(:)),Rithm(:),num2cell(iter(:)),num2cell(flag(:))...
                  cellstr(char(DateDepo(:))),cellstr(char(DateArro(:))),num2cell(TOF1o(:)),num2cell(TOF2o(:)),...
                  num2cell(totTOFo(:)), num2cell(totTOFo(:)/30), num2cell(DV1o(:)),num2cell(DV2o(:))];
        tnames2 = ["Name","Spkid","Opt.Iter","Cost-Fcn","Constraints","Algorithm","No.Iterations","Exit-Flag",...
                  "Dep Date","Arr Date","TOF1 (days)","TOF2 (days)","TOF total (days)","TOF total (months)",...
                  "Dep.Vel (km/s)","Rel.Vel (km/s)"];
    case 2
        %Exporting Original Results
        OriDat = [cellstr(AstName(1,:)),num2cell(spkid(1,:)),empty,empty,empty,empty,empty,empty,cellstr(char(DateDepi)),cellstr(char(DateArri)),num2cell(TOF1i(:)),num2cell(TOF2i(:)),...
                  num2cell(TOF3i(:)),num2cell(totTOFi(:)),num2cell(totTOFi(:)/30),num2cell(DV1i(:)), num2cell(DV2i(:)),num2cell(totDVi(:)),num2cell(DV3i(:))];
        tnames1 = ["Name","Spkid","Opt.Iter","Cost-Fcn","Constraints","Algorithm","No.Iterations","Exit-Flag","Dep Date","Arr Date","TOF1 (days)","TOF2 (days)",...
                   "TOF3 (days)", "TOF total (days)","TOF total (months)","Dep.Vel 1 (km/s)","Dep.Vel 2 (km/s)","Tot.Dep.Vel (km/s)","Rel.Vel (km/s)"];
        %Exporting Optimsation Data
        OptDat = [AstName(:),spkid(:),OptCount(:),coststr(:)',num2cell(scenario(:)'),Rithm(:)',num2cell(iter(:)'),num2cell(flag(:)')...
                  cellstr(char(DateDepo(:)')),cellstr(char(DateArro(:)')),num2cell(TOF1o(:)'),num2cell(TOF2o(:)'),num2cell(TOF3o(:)')...
                  num2cell(totTOFo(:)'), num2cell(totTOFo(:)'/30), num2cell(DV1o(:)'),num2cell(DV2o(:)'),num2cell(totDVi(:)'),num2cell(DV3i(:)')];
        tnames2 = ["Name","Spkid","Opt.Iter","Cost-Fcn","Constraints","Algorithm","No.Iterations","Exit-Flag",...
                  "Dep Date","Arr Date","TOF1 (days)","TOF2 (days)","TOF3 (days)","TOF total (days)","TOF total (months)",...
                  "Dep.Vel 1 (km/s)","Dep.Vel 2 (km/s)","Tot.Dep.Vel (km/s)","Rel.Vel (km/s)"];
    case 4
        %Exporting Original Results
        MM = table2cell(AlignedTable(month,1));
        OriDat = [MM(1,:),cellstr(AstName(1,:)),num2cell(spkid(1,:)),empty,empty,empty,empty,empty,empty,cellstr(char(DateDepi)),cellstr(char(DateArri)),num2cell(TOF1i(:)),num2cell(TOF2i(:)),...
                  num2cell(TOF3i(:)),num2cell(TOF4i(:)),num2cell(totTOFi(:)),num2cell(totTOFi(:)/30),num2cell(DV1i(:)), num2cell(DV2i(:)),num2cell(totDVi(:)),num2cell(DV3i(:))];
        tnames1 = ["Month","Name","Spkid","Opt.Iter","Cost-Fcn","Constraints","Algorithm","No.Iterations","Exit-Flag","Dep Date","Arr Date","TOF1 (days)","TOF2 (days)",...
                   "TOF3 (days)","TOF4 (days)", "TOF total (days)","TOF total (months)","Dep.Vel 1 (km/s)","Dep.Vel 2 (km/s)","Tot.Dep.Vel (km/s)","Rel.Vel (km/s)"];
        %Exporting Optimsation Data
        MM = repmat(table2cell(AlignedTable(month,1)),length(OptCount),1);
        OptDat = [MM(:),AstName(:),spkid(:),OptCount(:),coststr(:),num2cell(scenario(:)),Rithm(:),num2cell(iter(:)),num2cell(flag(:))...
                  cellstr(char(DateDepo(:))),cellstr(char(DateArro(:))),num2cell(TOF1o(:)),num2cell(TOF2o(:)),num2cell(TOF3o(:)),num2cell(TOF4o(:))...
                  num2cell(totTOFo(:)), num2cell(totTOFo(:)/30), num2cell(DV1o(:)),num2cell(DV2o(:)),num2cell(totDVi(:)),num2cell(DV3i(:))];
        tnames2 = ["Month","Name","Spkid","Opt.Iter","Cost-Fcn","Constraints","Algorithm","No.Iterations","Exit-Flag",...
                  "Dep Date","Arr Date","TOF1 (days)","TOF2 (days)","TOF3 (days)","TOF4 (days)","TOF total (days)","TOF total (months)",...
                  "Dep.Vel 1 (km/s)","Dep.Vel 2 (km/s)","Tot.Dep.Vel (km/s)","Rel.Vel (km/s)"];
end
%Displaying original results
OriTable = array2table(OriDat,"VariableNames",tnames1);
disp(OriTable)
%Displaying optimised results
OptTable = array2table(OptDat,"VariableNames",tnames2);
disp(OptTable)
switch Export
    case 1
        switch Design
            case 1
                writetable(OriTable,'OptimisationResults.xlsx','sheet','MS1','WriteMode','Append','WriteRowNames',true,'WriteVariableNames',true);
                writetable(OptTable,'OptimisationResults.xlsx','sheet','MS1','WriteMode','Append','WriteRowNames',true);
            case 2
                writetable(OriTable,'OptimisationResults.xlsx','sheet','MS2','WriteMode','Append','WriteRowNames',true,'WriteVariableNames',true);
                writetable(OptTable,'OptimisationResults.xlsx','sheet','MS2','WriteMode','Append','WriteRowNames',true);
            case 4
                writetable(OriTable,'OptimisationResults.xlsx','sheet','MS4','WriteMode','Append','WriteRowNames',true,'WriteVariableNames',true);
                writetable(OptTable,'OptimisationResults.xlsx','sheet','MS4','WriteMode','Append','WriteRowNames',true);
        end
end

%% Storing new transfer point 
options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16); %Setting tolerances (varied for situation)
[TOF1opt,Y1opt] = ode45(@(t,X)CR3BPEOM(mubar,X,t),[0 Xopt(7)],Xopt(1:6),options); %Integrating system of 1st order ODEs
OptimalTrajectory = sprintf('D%dOptimalLeg1%s.mat',Design,string(AstName));
save(OptimalTrajectory,'Y1opt','TOF1opt')

%% Clearing 'ans' in Command Window
if nargout<1 
    clear Xopt
end
end

