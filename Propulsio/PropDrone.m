clear
%% Torque Requirements for a Rectangular Propeller
% Propeller physical caratheristics
radius = 0.5; %[m] Distance from the hub to the tip
chord = 0.12; %[m] Assumed constant chord
pitch = [14 7]; %[º] Angle between the airfoil's chord and the hub's plane


% Airfoil selected - S1223-IL
Re = [50000 100000 200000 500000 1000000];

% Section calculations
g = 9.81; %[m/s^2] Gravity Acceleration
rho = 1.225; %[kg/m^3] Air Density
mu = 1.8e-5; %[Ns/m] Dynamic Viscosity
rpm = 3500; %[rmp] Propeller Turn-speed
omega = rpm*2*pi/60; %[rad/s]
elements = 100; %Number of domain elements
S = chord*(radius/elements); %[m^2] Element Surface


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

%Profile Alphas
alpha = linspace(pitch(1),pitch(2), elements);
%%
%Aerodynamic Variables
Total_Lift = 0;
Total_Drag = 0;
Total_Torque = 0;

 Cl = zeros(elements,1);
 Cd = zeros(elements,1);
 Re_X= zeros(elements,1);

for i = 1:elements
    x = radius*i/elements; % Central position of the element
    Re_X(i,1) =(rho*omega*x*chord)/mu;
   
    
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
    
    
    
    

    Total_Lift = Total_Lift + 0.5*rho*S*(omega*x)^2*Cl(i,1);
    Total_Drag = Total_Drag + 0.5*rho*S*(omega*x)^2*Cd(i,1);
    Total_Torque = Total_Torque + 0.5*rho*S*(omega*x)^2*Cd(i,1)*x;
    
end


% Double Bladed Propeller
Total_Lift = 2*Total_Lift;
Total_Drag = 2*Total_Drag;
Total_Torque = 2*Total_Torque;
% Units Adaptation
Thrust = Total_Lift/g;
Power = Total_Torque*omega;
%He fet aquesta tonteria per intentar de donar numeros ja de Torque, Power etc
%Cal dir que no son massa grans els valors que surten. Com a molt 5Kw