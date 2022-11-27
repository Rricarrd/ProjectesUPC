clc
clear

%% ################## Drone Axial Flight ##########################
%% Initial variables
%Some constants
g = 9.81; %[]
pi = 3.141592;

% Atmospherical conditions
z = 0; [rho,~,T,a,visco_din] = AtmosphereConditions(z);

%% Rotor parameters
%Rotor radius
R = 1; %[m]
A = pi*R^2; %[m^2]

%Blade number
Nb = 2; %[]

% Rotor solidity
sigma = 0.1; %[Relationship of total blade area vs rotor area  ]

% Twist parameters (Lineal twist)
theta_tw_deg = -5; %[deg/m]
theta_tw = DegToRad(theta_tw_deg); %[rad/m]

%Airfoil parameters
CL_alpha_deg = 0.11; %[1/deg]
% CL_alpha = 1/DegToRad(1/CL_alpha_deg); %[1/rad]
CL_alpha = 2*pi;
Cd0 = 0.01;

%Thrust coefficients
MTOW = 520; %[kg]
W = MTOW*g; %[N]
RPMs = 2500;
omega = RPMs*2*pi/60;
Vtip = omega*R;

% CTreq = (W/8)/(rho*A*Vtip^2);
CTreq = 0.008;
%% Numerical parameters
error = 0.01;
elements = 50;
deltar = R/elements;
r = linspace(deltar/2,R,elements);

%Initial definitions
deltaCT = zeros(elements-1,1);
deltaCPi = zeros(elements-1,1);
deltaCP0 = zeros(elements-1,1);
theta_0 = ((6*CTreq)/(sigma*CL_alpha))-((3/4)*theta_tw) + ((3/2)*sqrt(CTreq/2));
convThetaOCT = false;

while convThetaOCT == false
%Some definitions
CT = 0;
CP = 0;    


    for j = 1:elements-1
    
    %Element pitch angle
    thetaj = theta_0 + r(j)*theta_tw;
        
    %Parameters
    Fk = 1;
    convF = false;
    
        % F Factor Iteration
        while convF == false
       
        %Inflow ratio
        lambdai = ((sigma*CL_alpha)/(16*Fk))*(-1+sqrt(1+((32*Fk*r(j)*thetaj)/(sigma*CL_alpha))));
        
        %Incidence angle
        phi = lambdai/r(j);
        
        %Prandtl Correction Parameter
        f = 0.5*Nb*((1-r(j))/(r(j)*phi));
        
        %Prandtl Correction Factor
        Fk_star = (2/pi)*(acos(exp(-f)));
        
            if abs(Fk - Fk_star) > error
            Fk = Fk_star;
            else
            convF = true;
            Fk = Fk_star;
            end   
        end
    
    %Coeficients
    deltaCT(j,1) = 4*Fk*lambdai^2*r(j);
    Cl = CL_alpha*(thetaj - lambdai/r(j));
    deltaCPi(j,1) = lambdai*deltaCT(j,1);
    deltaCP0(j,1) = 0.5*sigma*Cd0*R^3;

    %Coeficients Summnation
    CT = CT + deltaCT(j,1)*r(j);

    end

theta_0_star = theta_0 + ((6*(CTreq-CT))/(sigma*CL_alpha)) + (3*sqrt(2)/4)*(sqrt(CTreq)-sqrt(CT));

if abs(CT-CTreq) > error || abs(theta_0-theta_0_star) > error
    theta_0 = theta_0_star;      
else
   convThetaOCT = true;
  
end


end

theta_0_deg = (theta_0/(2*pi))*360;

%Power coefficient

for j = 1:elements-1    
CP = CP + deltaCPi(j,1) + deltaCP0(j,1)*r(j);
end

P = CP*rho*A*Vtip^3; %[W]
P_kW = P/1000;
