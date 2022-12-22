clear
clc

%% Mass and volume calculation for the battery system
Pm = 12330; % potencia 1 motor [W]
N = 8; % number of motors
Pmot=Pm*N; % total power for the motors
Pav = 1000; %[W] ??


Ft_min= 1:1:60; %Flight time [min]
Ft = Ft_min/60; %Flight time [h]

Emot = Pmot*Ft; %[Wh]
Eav = Pav*Ft; %[Wh]
E = (Emot+Eav)/1000; %[kWh]

%% Amprius cells data
P_W = 450 ; % power/weight ratio[Wh/kg]

P_Vol = 1250e3; % energy density [Wh/L]


%% Calculations
%Masses
Mmot = Emot/P_W; % mass of the battery system for the motors [kg]
Mav = Eav/P_W; % mass of the battery system for the avionics [kg]

%Volumes
Vmot = Emot/P_Vol; % volume of the battery system for the motors [m^3]
Vav = Eav/P_Vol;  % volume of the battery system for the avionics [m^3]
Vmot_l = Vmot*1000; % volume of the battery system for the motors [l]
Vav_l = Vav*1000; % volume of the battery system for the avionics [l]

%Cubical dimensions approx
cube_side = nthroot(Vmot, 3)*100; % side of a cube that contains the batteries [cm]

%% Plot
plot(Ft_min,Mmot+Mav);
xlabel('Flight time (minutes)');
ylabel('Battery weight (kg)');
title('Battery weight versus flight time diagram ')

coefs = polyfit(Ft_min,Mmot+Mav,1);
pendent = coefs(1);


%% Flight time / payload graph
MTOW = 518.75; %[kg]
PL = 1:1:280; %[kg]
BAT = 140; %[kg]
OEW_nobat = 80; %[kg]




% % Size calculation
% zmax = 0.2; % allowed depth for the motor battery system [m]
% area_mot = Vmot/zmax; % alllowed area for the motor battery system [m^2]
% x_mot = sqrt(area_mot); 
% x_mot_div=x/8; % size battery system cases for the motors
% 
% zav = 0.1; % allowed depth for the avionics battery system [m]
% area_av = ; % alllowed area for the avionics battery system [m^2]
% x_av = sqrt(area_av); % size battery system case for the avionics





% syms xmot ymot xav yav
% 
% eq1 = Vmot == zmax*xmot*ymot;
% eq2 = 2*xmot + 2*ymot == P; 
% %eq3 = Vav == zav*xav*yav;
% %eq4 = 2*xav + 2*yav == P_av; 
% [sol1, sol2, sol3, sol4] = vpasolve([eq1 eq2 eq3 eq4], [xmot ymot xav yav]) ;
% xmot_sol = sol1;
% ymot_sol = sol2;
% %xav_sol = sol3;
% %yav_sol = sol4;
















