clear
clc
% Torque Requirements for a Rectangular Propeller
% Propeller physical caratheristics
radius = 0.5; %[m] Distance from the hub to the tip ########### FIXAT
chord = [0.12 0.08]; %[m] Assumed constant chord ########### ES POT VARIAR
pitch = [14 7]; %[º] Angle between the airfoil's chord and the hub's plane ########### ES POT VARIAR
n_blades = 2;  %########### ES POT VARIAR

% Airfoil selected - S1223-IL
Re = [50000 100000 200000 500000 1000000];

% Section calculations
g = 9.81; %[m/s^2] Gravity Acceleration
rho = 1.225; %[kg/m^3] Air Density
mu = 1.8e-5; %[Ns/m] Dynamic Viscosity
rpm = 3000; %[rmp] Propeller Turn-speed   ########### ES POT VARIAR
omega = rpm*2*pi/60; %[rad/s]
elements = 100; %Number of domain elements
pi = 3.141592;
oswald = 0.85;

%Import tables
[Re5e4,Re1e5,Re2e5,Re5e5,Re1e6] = importTables()

%