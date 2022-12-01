clc
clear
%% Some calculations
%Global Power 
P = 100; %[kW]
V = 70; %[V]
I = P/V *1000; %[A]

%Motor Power 
Pmot = 15; %[kW]
Vmot = 70; %[V]
Imot = Pmot/Vmot *1000; %[A]

%Material properties
Al_cond = 3.5e7; %[S/m]
Al_res = 1/Al_cond; %[ohm*m]
Al_dens = 2700; % [kg/m3]

Cu_cond = 5.98e7; %[S/m]
Cu_res = 1/Cu_cond; %[ohm*m]
Cu_dens = 8960; % [kg/m3]

%Copper mass to aluminium mass (Same R and L, different A)
mCu_mAl = (Cu_res*Cu_dens)/(Al_res*Al_dens);




%Main buses dimensioning
Ploses = 30; %[W]
Ibar = I; %Considering the whole bar doesn't need to support all the Amps and the distances are smaller
Rreq = Ploses/(Ibar^2); %[Ohm]
L = 0.2; %[m]
S = (Al_res*L)/Rreq; %[m^2]
Scm = S*100^2; %[cm^2]
c = sqrt(Scm); %[cm]
mbar = Al_dens*S*L*4;
mbuses = mbar*2;


%Motor cables dimensioning
Ploses_cables = 5; %[W]
Icables = Imot; %Considering the whole bar doesn't need to support all the Amps and the distances are smaller
Rreq_cables = Ploses_cables/(Icables^2); %[Ohm]
Lcables = 0.5; %[m]
Scables = (Al_res*L)/Rreq_cables; %[m^2]
Scables_cm = Scables*100^2; %[cm^2]
rcm = sqrt(Scables_cm/pi); %[cm]
mcable = Al_dens*Scables*Lcables;
mcables = mcable*16;
