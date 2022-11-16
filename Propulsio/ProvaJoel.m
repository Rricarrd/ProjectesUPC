%% Torque Requirements for a Rectangular Propeller
% Propeller physical caratheristics
radius = 0.5; %[m] Distance from the hub to the tip
chord = 0.15; %[m] Assumed constant chord
pitch = 7.5; %[º] Angle between the airfoil's chord and the hub's plane
% Airfoil selected - NACA 0012
Re = [50000 100000 200000 500000 1000000];
% cl = [0.65 0.70 0.75 0.80 0.90]; % Original
cl = [2 2 2 2 2]; % If cl was constant from s-1223-il
cd = [0.6 0.4 0.3 0.2 0.1];
% Section calculations
g = 9.81; %[m/s^2] Gravity Acceleration
rho = 1.225; %[kg/m^3] Air Density
mu = 1.8e-5; %[Ns/m] Dynamic Viscosity
rpm = 2500; %[rmp] Propeller Turn-speed
omega = rpm*2*pi/60; %[rad/s]
elements = 100; %Number of domain elements
S = chord*(radius/elements); %[m^2] Element Surface
Lift = zeros(elements,1);
 Drag = zeros(elements,1);
 Torque = zeros(elements,1);
for i = 1:elements
    x = radius*i/elements; % Central position of the element
    Re_x =(rho*omega*x*chord)/mu;
    Cl = 0;
    Cd = 0;
    for j =2:numel(Re)
        if Re_x<Re(j)
        Cl = cl(j);
        Cd = cd(j);
        break;
        elseif Re_x<1000000000
        Cl = 2;
        Cd = 0.1;
        break
        end
    end
     
    Lift(i,1) = 0.5*rho*S*(omega*x)^2*Cl;
    Drag(i,1) = 0.5*rho*S*(omega*x)^2*Cd;
    Torque(i,1) = 0.5*rho*S*(omega*x)^2*Cd*x;

    Total_Lift = Total_Lift + Lift(i,1) ;
    Total_Drag = Total_Drag +  Drag(i,1) ;
    Total_Torque = Total_Torque +  Torque(i,1);
end
% Double Bladed Propeller
Total_Lift = 2*Total_Lift;
Total_Drag = 2*Total_Drag;
Total_Torque = 2*Total_Torque;
% Units Adaptation
Thrust = Total_Lift/g;
Power = Total_Torque*omega;