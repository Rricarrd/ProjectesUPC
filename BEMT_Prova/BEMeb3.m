function [cl,cd,phi,aoa,a,ap]=BEMeb3
%% Preallocate space for variables within loop
%Determine size of test vectors/arrays
ns=length([5 8 11 14 17 21]);
RNodes=[0.2 1.0 2.0 3.0 4.0 5.0];
Chord=[0.7 0.71 0.44 0.3 0.23 0.19];
twst=pi/180*[29.0 15.7 5.1 0.9 -1.3 -2.6];

rhub=0.2;
rtip=5.0;
a0=zeros(1,ns); %Old (previous iteration) axial induction factor
ap0=zeros(1,ns); %Old (previous iteration) tangential induction factor
phi=zeros(1,ns); %Local inflow angle
aoa=zeros(1,ns); %Local angle of attack
cl=zeros(1,ns); %Local lift coefficient
cd=zeros(1,ns); %Local drag coefficient
ct=zeros(1,ns); %Local thrust coefficient
ftip=zeros(1,ns); %Tip loss factor
fhub=zeros(1,ns); %Hub loss factor
f=zeros(1,ns); %Total loss corection factor
fiter=zeros(1,ns); %Converence flag for gridpoints ('1' if converged, '9999' if not)

%% Define convergence criteria
tol=1e-6; %Convergence tolerance
da=ones(1,ns); %Set initial value for axial induction factor residual equal to 1
dap=ones(1,ns); %Set initial value for tangential induction factor residual equal to 1
ncv=find(da>tol | dap>tol); %Identify all nonconverged points (all initially)
miter=5000; %Maximum number of allowable iterations
wt=0.5; %Weighting factor on corrections to balance speed with stability (faster as you 
        %approach 1, but less stable)

%% Compute relevant velocity/angle components
Uinf=4.5;
Om=11.2;
ptch=0.0;

%% Compute initial guesses of key variables
ptch=ptch';
sigmap=3*Chord./(2.*pi.*RNodes); %Local solidity
lambdar=Om*RNodes/Uinf; %Local speed ratio

% Initial values for axial and tangential induction factors
a=real(0.25*(2+pi*lambdar.*sigmap-sqrt(4-4*pi*lambdar.*sigmap+pi*lambdar.^2.*sigmap.*(8*(twst+ptch)+pi*sigmap))));
size(a)
size(a0)
ap=zeros(size(a));
a=zeros(size(a));
%% Primary loop for BEM
for j=1:3000
    % Save previous values of axial and tangential induction factors
    a0(ncv)=a(ncv);
    ap0(ncv)=ap(ncv);
    
    % Compute inflow angle and angle of attack
    a
    ap
    phi(ncv)=atan2(Uinf*(1-a(ncv)),Om*RNodes(ncv).*(1+ap(ncv)));
    aoa(ncv)=(phi(ncv)-(twst(ncv)+ptch));
    
    % Interpolate over airfoil database for lift and drag coefficients
for ii=1:ns
    if (aoa(ii)<0.2443)
    cl(ii)=26.95*aoa(ii)^3-16.67*aoa(ii)^2 + 7.9*aoa(ii)-0.02195;
    end
    if (aoa(ii)>0.2443)
    cl(ii)=(-0.1275*aoa(ii)^2 +0.2109*aoa(ii)-0.0388)/(aoa(ii)^3-1.721*aoa(ii)^2+0.9863*aoa(ii)-0.1495);
    end
   
    if (aoa(ii)<0.1745)
    cd(ii)=0.0;
    end
    if (aoa(ii)>0.1745)
    cd(ii)=0.0;
    end   
   
end
  cl
  cd(ncv)=0.0;
    % Compute elemental thrust coefficient
    ct(ncv)=sigmap(ncv).*(1-a(ncv)).^2.*(cl(ncv).*cos(phi(ncv))+cd(ncv).*sin(phi(ncv)))./ ...
        sin(phi(ncv)).^2;
    
    % Compute loss correction factor due to tip and hub losses
    f(ncv)=1.0;
    % Compute axial induction factor using conventional BEM theory

    a(ncv)=((1+4.*f(ncv).*sin(phi(ncv)).^2./(sigmap(ncv).*(cl(ncv).*cos(phi(ncv))+cd(ncv).*sin(phi(ncv))))).^-1);
      % Identify highly loaded gridpoints (requires use of modified Glauert correction for 
    % axial induction factor)
    
        % Compute tangential induction factor
    ap(ncv)=sigmap(ncv).*cl(ncv)./(4.0*(Om*RNodes(ncv)/Uinf).*sin(phi(ncv))).*(1.0-a(ncv));
       
    % Compute residuals
    da(ncv)=abs(a0(ncv)-a(ncv));
    dap(ncv)=abs(ap0(ncv)-ap(ncv));
    
    % Apply corrective weighting for convergence stability
    if wt>0
        a(ncv)=a0(ncv)+wt.*(a(ncv)-a0(ncv));
        ap(ncv)=ap0(ncv)+wt.*(ap(ncv)-ap0(ncv));
    end

    % Clear all gridpoint flags in preparation for next loop
    clear ncv ncvf ncvcl ida idap
    
    % Identify nonconverged gridpoints
    ncv=find(da>tol | dap>tol);
    
    % If all points meet convergence criteria, break loop
    if isempty(ncv)
        break
    end
    
    % If maximum allowable iterations has been reached, flag nonconverged gridpoints 
    % with '9999'
    if j==miter
        fiter(ncv)=9999;
    else
        fiter(ncv)=j;
    end
end

lambdar
ncv
plot(RNodes,a)
figure
plot(RNodes,ap)
ef=lambdar.^3.*(1-a).*ap.*(1-cd./cl.*cot(phi)).*f;
dela=(lambdar(2:ns)-lambdar(1:ns-1))/2;
cipi=8/(lambdar(ns))^2*sum(dela.*(ef(2:ns)+ef(1:ns-1)))