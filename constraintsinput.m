function [c,ceq] = constraintsinput(Design,X,scenario,Event,Yend,Y1new,Y2new,RV_Ast,ropt1,vopt1,ropt2,vopt2,ropt3,mubar,TU,LU,mu)

c = []; %Default if no other c's selected

switch Design
    case 1
        %Fundamentals
        ceq(1:3) = Yend(1:3) - RV_Ast(1:3); %Fundamental constraint - Asteroid Reached 
        ceq(4:6) = Yend(4:6) + X(12:14) - RV_Ast(4:6); %Relative velocity evens out
        %Events
        r1 = [ropt1(1) + mubar ; ropt1(2:3)]; %Position vector with respect to the Sun   
    case 2
        %Fundamentals
        ceq(1:3) = Yend(1:3) - RV_Ast(1:3)/LU; %Fundamental constraint - Asteroid Reached 
        ceq(4:6) = Yend(4:6) + X(16:18) - RV_Ast(4:6)/LU*TU; %Relative velocity evens out 
        %Events
        r1 = [ropt1(1)-(1-mubar) ; ropt1(2:3)]; %Position vector with respect to the Earth
        %Defining SOI boundary parameters 
        Msun = 1.989e30; Mearth =5.972e24;
        SOI = 0.9431*(Mearth/Msun)^(2/5); %Normalised SOI radius
        r2 = [ropt2(1)-(1-mubar) ; ropt2(2:3)]; %Position vector with respect to the Earth
        c(1) = SOI - norm(r2); %S/C is outside Earth SOI - also fundamental
        %Defining eccentricity after burn 1 at Perigee for hyperbolic constraint
        ropt1burn = Y1new(1:3); vopt1burn = Y1new(4:6);
        r1burn = [ropt1burn(1)-(1-mubar) ; ropt1burn(2:3)]; %Position vector wrt Earth
        h = cross(r1burn,vopt1burn); %Angular Momentum
        e_vector = (1/mu*cross(vopt1burn,h))-r1burn/norm(r1burn); % Eccentricity vector
        e1 = norm(e_vector); %Eccentricity at Perigee after burn
    case 4
        %Fundamentals
        ceq(1:3) = Yend(1:3) - RV_Ast(1:3); %Fundamental constraint - Asteroid Reached 
        ceq(4:6) = Yend(4:6) + X(17:19) - RV_Ast(4:6); %Relative velocity evens out
        %Events
        r2 = [ropt2(1)-(1-mubar) ; ropt2(2:3)]; %Position vector with respect to the Earth
        %Defining SOI boundary parameters 
        Msun = 1.989e30; Mearth =5.972e24;
        SOI = 0.9431*(Mearth/Msun)^(2/5); %Normalised SOI radius
        r3 = [ropt3(1)-(1-mubar) ; ropt3(2:3)]; %Position vector with respect to the Earth
        c(1) = SOI - norm(r3); %S/C is outside Earth SOI - also fundamental
        %Defining eccentricity after burn 1 at Perigee for hyperbolic constraint
        ropt2burn = Y2new(1:3); vopt2burn = Y2new(4:6);
        r2burn = [ropt2burn(1)-(1-mubar) ; ropt2burn(2:3)]; %Position vector wrt Earth
        h = cross(r2burn,vopt2burn); %Angular Momentum
        e_vector = (1/mu*cross(vopt2burn,h))-r2burn/norm(r2burn); % Eccentricity vector
        e1 = norm(e_vector); %Eccentricity at Perigee after burn
end

switch Design
    case 1
        if Event == 1
            ceq(7) = dot(r1,vopt1); %Entry into Aphelion
        end
        %Non-linear inequality constraints
        if scenario == 12 
            c(1) = norm(X(8:10)) - 0.15*TU/LU; %Searching for departure velocity compliance 
            c(2) = norm(X(12:14)) - 5.0*TU/LU; %Searching for arrival (relative) velocity compliance 
        elseif scenario == 1
            c(1) = norm(X(8:10)) - 0.15*TU/LU; %Searching for departure velocity compliance 
        elseif scenario == 2
            c(1) = norm(X(12:14)) - 5.0*TU/LU; %Searching for arrival (relative) velocity compliance 
        end

    case 2
        if Event == 1
            ceq(7) = dot(r1,vopt1); %Entry into Perigee
        end
    %Non-linear inequality constraints
    if scenario == 123
        c(2) = 1 - e1; %Flyby at Perigee after burn is unbounded 
        c(3) = (norm(X(8:10)) + norm(X(12:14)))- 0.15*TU/LU; %Searching for departure velocity compliance 
        c(4) = norm(X(16:18)) - 5.0*TU/LU; %Searching for arrival (relative) velocity compliance 
    elseif scenario == 12
        c(2) = 1 - e1; %Flyby at Perigee after burn is unbounded 
        c(3) = (norm(X(8:10)) + norm(X(12:14)))- 0.15*TU/LU; %Searching for departure velocity compliance
    elseif scenario == 13
        c(2) = 1 - e1; %Flyby at Perigee after burn is unbounded 
        c(3) = norm(X(16:18)) - 5.0*TU/LU; %Searching for arrival (relative) velocity compliance 
    elseif scenario == 23
        c(2) = (norm(X(8:10)) + norm(X(12:14)))- 0.15*TU/LU; %Searching for departure velocity compliance 
        c(3) = norm(X(16:18)) - 5.0*TU/LU; %Searching for arrival (relative) velocity compliance 
    elseif scenario == 1
        c(2) = 1 - e1; %Flyby at Perigee after burn is unbounded 
    elseif scenario == 2
        c(2) = (norm(X(8:10)) + norm(X(12:14)))- 0.15*TU/LU; %Searching for departure velocity compliance %(4)
    elseif scenario == 3
        c(2) = norm(X(16:18)) - 5.0*TU/LU; %Searching for arrival (relative) velocity compliance %(5)
    end

    case 4
        if Event == 1
            ceq(7) = dot(r2,vopt2); %Entry into Perigee
        end
    %Non-linear inequality constraints
    if scenario == 123
        c(2) = 1 - e1; %Flyby at Perigee after burn is unbounded 
        c(3) = (norm(X(9:11)) + norm(X(13:15)))- 0.15*TU/LU; %Searching for departure velocity compliance 
        c(4) = norm(X(17:19)) - 5.0*TU/LU; %Searching for arrival (relative) velocity compliance 
    elseif scenario == 12
        c(2) = 1 - e1; %Flyby at Perigee after burn is unbounded 
        c(3) = (norm(X(9:11)) + norm(X(13:15)))- 0.15*TU/LU; %Searching for departure velocity compliance 
    elseif scenario == 13
        c(2) = 1 - e1; %Flyby at Perigee after burn is unbounded 
        c(4) = norm(X(17:19)) - 5.0*TU/LU; %Searching for arrival (relative) velocity compliance
    elseif scenario == 23
        c(3) = (norm(X(9:11)) + norm(X(13:15)))- 0.15*TU/LU; %Searching for departure velocity compliance  
        c(4) = norm(X(17:19)) - 5.0*TU/LU; %Searching for arrival (relative) velocity compliance 
    elseif scenario == 1
        c(2) = 1 - e1; %Flyby at Perigee after burn is unbounded 
    elseif scenario == 2
        c(3) = (norm(X(9:11)) + norm(X(13:15)))- 0.15*TU/LU; %Searching for departure velocity compliance 
    elseif scenario == 3
        c(4) = norm(X(17:19)) - 5.0*TU/LU; %Searching for arrival (relative) velocity compliance
    end
end