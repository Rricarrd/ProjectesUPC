clear
clc








%% ########################################################
%% Thrust and Power Calculation
% Propeller physical caratheristics
radius = 0.5; %[m] Distance from the hub to the tip ########### FIXAT
chord = [0.08 0.025]; %[m] Assumed constant chord ########### ES POT VARIAR
pitch = [14 3]; %[º] Angle between the airfoil's chord and the hub's plane ########### ES POT VARIAR
n_blades = 2;  %########### ES POT VARIAR
rpm_min = 1000;
rpm_max = 5500;

% Airfoil selected - S1223-IL
Re = [50000 100000 200000 500000 1000000];

% Section calculations
g = 9.81; %[m/s^2] Gravity Acceleration
rho = 1.225; %[kg/m^3] Air Density
mu = 1.8e-5; %[Ns/m] Dynamic Viscosity


elements = 100; %Number of domain elements
pi = 3.141592;
oswald = 0.85;
motor_efficiency = 0.9;
ESC_efficiency = 0.9;

%Airfoil Data S1223-IL
Re5e4_tab = readtable('xf-s1223-il-50000-n5.csv');
Re1e5_tab = readtable('xf-s1223-il-100000-n5.csv');
Re2e5_tab = readtable('xf-s1223-il-200000-n5.csv');
Re5e5_tab = readtable('xf-s1223-il-500000-n5.csv');
Re1e6_tab = readtable('xf-s1223-il-1000000-n5.csv');

Re5e4 = table2array(Re5e4_tab);
Re1e5 = table2array(Re1e5_tab);
Re2e5 = table2array(Re2e5_tab);
Re5e5 = table2array(Re5e5_tab);
Re1e6 = table2array(Re1e6_tab);

%Profile Chords
local_chord = linspace(chord(1),chord(2), elements);

%Profile Alphas
alpha = linspace(pitch(1),pitch(2), elements);

%Aerodynamic Variables
 Lift = zeros(elements,1);
 Drag = zeros(elements,1);
 Torque = zeros(elements,1);
 Efficiency_wing = zeros(elements,1);
 Efficiency_airfoil = zeros(elements,1);

 CD = zeros(elements,1);
 CL = zeros(elements,1);
 Cl = zeros(elements,1);
 Cd = zeros(elements,1);
 Re_X= zeros(elements,1);
 
Total_Lift = zeros(elements,1); %[N]
Total_Drag = zeros(elements,1); %[N]
Total_Torque = zeros(elements,1); %[Nm]

THRUST= zeros(elements,1); %[kgf]
MECHANICAL_POWER = zeros(elements,1); %[kW]
ELECTRICAL_POWER = zeros(elements,1);
TP_RATIO = zeros(elements,1);%[kgf/kW]


rpm = linspace(rpm_min,rpm_max,elements);


for k = 1:elements

omega = (rpm(1,k)*2*pi)/60; %[rad/s]

for i = 1:elements
    
    %Local Reynolds
    x = radius*i/elements; % Central position of the element
    Re_X(i,1) =(rho*omega*x*local_chord(i))/mu;
   
    
    if Re_X(i,1) <= Re(1) %%Re<50000 ######################
       
     for j = 1:size(Re5e4)-1
      if Re5e4(j+1,1)>alpha(i) && Re5e4(j,1)<=alpha(i)
      Cl2 = Re5e4(j,2);   
      Cd2 = Re5e4(j,3);     
      end 
     end
     
     Cl1 = 0;
     Cd1 = 0;
     
     Cl(i,1) = Interpol(0,Cl1,Re(1),Cl2,Re_X(i,1));
     Cd(i,1) = Interpol(0,Cd1,Re(1),Cd2,Re_X(i,1));
     
    elseif  Re_X(i,1)>Re(1) && Re_X(i,1)<=Re(2) %% 50000<Re<100000  ######################
       
     for j = 1:size(Re5e4)-1
      if Re5e4(j+1,1)>alpha(i) && Re5e4(j,1)<=alpha(i)
      Cl1 = Re5e4(j,2);   
      Cd1 = Re5e4(j,3);     
      end 
     end
     
     for j = 1:size(Re1e5)-1
      if Re1e5(j+1,1)>alpha(i) && Re1e5(j,1)<=alpha(i)
      Cl2 = Re1e5(j,2);   
      Cd2 = Re1e5(j,3);     
      end 
     end
        
    Cl(i,1) = Interpol(Re(1),Cl1,Re(2),Cl2,Re_X(i,1));
    Cd(i,1) = Interpol(Re(1),Cd1,Re(2),Cd2,Re_X(i,1)); 
     
    elseif  Re_X(i,1)>Re(2) && Re_X(i,1)<=Re(3)   %% 100000<Re<200000  ######################
      
     for j = 1:size(Re1e5)-1
      if Re1e5(j+1,1)>alpha(i) && Re1e5(j,1)<=alpha(i)
      Cl1 = Re1e5(j,2);   
      Cd1 = Re1e5(j,3);     
      end 
     end
     
     for j = 1:size(Re2e5)-1
      if Re2e5(j+1,1)>alpha(i) && Re2e5(j,1)<=alpha(i)
      Cl2 = Re2e5(j,2);   
      Cd2 = Re2e5(j,3);     
      end 
     end
    
    Cl(i,1) = Interpol(Re(2),Cl1,Re(3),Cl2,Re_X(i,1));
    Cd(i,1) = Interpol(Re(2),Cd1,Re(3),Cd2,Re_X(i,1)); 
        
    elseif  Re_X(i,1)>Re(3) && Re_X(i,1)<=Re(4)   %% 200000<Re<500000  ######################
        
     for j = 1:size(Re2e5)-1
      if Re2e5(j+1,1)>alpha(i) && Re2e5(j,1)<=alpha(i)
      Cl1 = Re2e5(j,2);   
      Cd1 = Re2e5(j,3);     
      end 
     end
     
     for j = 1:size(Re5e5)-1  
      if Re5e5(j+1,1)>alpha(i) && Re5e5(j,1)<=alpha(i)
      Cl2 = Re5e5(j,2);   
      Cd2 = Re5e5(j,3);    
      end 
     end
     
    Cl(i,1) = Interpol(Re(3),Cl1,Re(4),Cl2,Re_X(i,1));
    Cd(i,1) = Interpol(Re(3),Cd1,Re(4),Cd2,Re_X(i,1)); 
     
      elseif  Re_X(i,1)>Re(4) && Re_X(i,1)<=Re(5)   %% 500000<Re<1000000  ######################
        
     for j = 1:size(Re5e5)-1  
      if Re5e5(j+1,1)>alpha(i) && Re5e5(j,1)<=alpha(i)
      Cl1 = Re5e5(j,2);   
      Cd1 = Re5e5(j,3);     
      end 
     end
     
     for j = 1:size(Re1e6)-1  
      if Re1e6(j+1,1)>alpha(i) && Re1e6(j,1)<=alpha(i)
      Cl2 = Re1e6(j,2);   
      Cd2 = Re1e6(j,3);     
      end 
     end
        
    Cl(i,1) = Interpol(Re(4),Cl1,Re(5),Cl2,Re_X(i,1));
    Cd(i,1) = Interpol(Re(4),Cd1,Re(5),Cd2,Re_X(i,1)); 
     
    elseif  Re_X(i,1)>Re(5)  %% Re<1000000  ######################
        
     for j = 1:size(Re1e6)-1  
      if Re1e6(j+1,1)>alpha(i) && Re1e6(j,1)<=alpha(i)
      Cl1 = Re1e6(j,2);   
      Cd1 = Re1e6(j,3);     
      end 
      
     Cl2 = 0;
     Cd2 = 0;
     end
     
      Cl(i,1) = Interpol(3e6,Cl2,Re(5),Cl1,Re_X(i,1));
      Cd(i,1) = Interpol(3e6,Cd2,Re(5),Cd1,Re_X(i,1)); 
      
    else 
    end
    
  
 
    %AR (induced drag assumed constant along the wing)
    if i<elements
    AR = (radius)^2/(0.5*(chord(1)+chord(2))*radius); %Whole Wing AR with trapezoidal area
    S = 0.5*(local_chord(i)+local_chord(i+1))*(radius/elements); %[m^2] Element Surface Trapezoidal
    else
    end
    
    %Coeficients de la pala completa
    CL(i,1) = Cl(i,1);
    CD(i,1) = Cd(i,1) + CL(i,1)^2/(AR*pi*oswald) ;
    
    
    %Lift i drag
    Lift(i,1) = 0.5*rho*S*(omega*x)^2*CL(i,1);
    Drag(i,1) = 0.5*rho*S*(omega*x)^2*CD(i,1);
    Torque(i,1) = 0.5*rho*S*(omega*x)^2*CD(i,1)*x;
    Efficiency_wing(i,1) = CL(i,1)/CD(i,1);
    Efficiency_airfoil(i,1) = Cl(i,1)/Cd(i,1);
    
    Total_Lift(k,1) = Total_Lift(k,1) + Lift(i,1) ;
    Total_Drag(k,1) = Total_Drag(k,1) +  Drag(i,1) ;
    Total_Torque(k,1) = Total_Torque(k,1) +  Torque(i,1);
    
end



% Double Bladed Propeller
Total_Lift(k,1) = n_blades*Total_Lift(k,1); %[N]
Total_Drag(k,1) = n_blades*Total_Drag(k,1); %[N]
Total_Torque(k,1) = n_blades*Total_Torque(k,1); %[Nm]

% Units Adaptation
THRUST(k,1) = (Total_Lift(k,1)/g); %[kgf]
MECHANICAL_POWER(k,1) = Total_Torque(k,1)*omega/1000; %[kW]
ELECTRICAL_POWER(k,1) = MECHANICAL_POWER(k,1)/(motor_efficiency*ESC_efficiency);
TP_RATIO(k,1) = THRUST(k,1)/ELECTRICAL_POWER(k,1);%[kgf/kW]


end


%% ######################### PL FT diagram
%% Data and polyfit

MTOW = 518.75; %[kg]
PL = 1:1:200; %[kg]
BAT = 140; %[kg]
OEW_nobat = 178; %[kg]

%Fixed weight
FixedW = OEW_nobat+BAT;


%Comparing polyfit to data
%Plot 
figure
plot(THRUST,ELECTRICAL_POWER);
xlabel('Thrust (kgf)');
ylabel('Power (kW))');
title('Thrust vs power');
hold on

%Polyfit
thrust_power = polyfit(THRUST,ELECTRICAL_POWER,3);
POWER_poly = polyval(thrust_power,THRUST);

plot(THRUST,POWER_poly);
hold off

%% Amprius cells data
P_W = 450 ; % power/weight ratio[Wh/kg]
P_Vol = 1250e3; % energy density [Wh/L]



%% Polyfit evaluation

Ft = zeros(1,length(PL));
P = zeros(1,length(PL));
for j = 1:length(PL)
    %Weight = Thrust
    Tm = (FixedW + PL(j))/8;
    
    
    %Motor power
    Pm = polyval(thrust_power,Tm);  % potencia 1 motor [W]  
    N = 8; % number of motors
    Pmot=Pm*N; % total power for the motors
    Pav = 1; %[W] ??
    P(j) = Pmot+Pav;
   
    
    %Battery Capacity
    FT_at_max_PL = 40/60; %[h]
    Pmax_PL = 12.33*8; %[kWh]
    Emot = Pmax_PL*FT_at_max_PL;
    Eav = Pav*FT_at_max_PL;
    E = Emot+Eav;
    

    Ft(j) = E/P(j); %[kWh]
    
end

figure
plot(PL,Ft*60)
xlabel('Payload Weight (kg)');
ylabel('Flight Time (min)');
title('Payload vs Flight Time')
