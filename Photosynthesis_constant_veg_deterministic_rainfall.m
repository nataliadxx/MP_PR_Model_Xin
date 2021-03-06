%--------------------------------------------------------------------------
%Phytoremediation_Modeling.m: models leaching of heavy metals and other
%contaminants from soils as well as plant-water uptake
%
%Authors: Xin Wang and Gabriel Katul
%
%Date: May 21, 2019
%
%References: 
%Manzoni et al. (2011) Proceedings of the Royal Society A, 467, 3188-3205
%
%--------------------------------------------------------------------------
%% The model
clear all
clc
%------------ Construct synthetic precipitation time series
freq=1;                       % Return frequency between days (1/d)
annual_precip=800;               % Annual precip. (mm/year)
N=365;                           % Number of days to simulate the process
dep=(annual_precip/365)/freq;
Pr=ones(1,N)*dep;                  %constant rainfall (mm/d)
tday=[0:1:N-1];                    % Time series
%----------- Plot synthetic precipitation time series
%Plot_Precip_Time_Series

Eff=[];
s=[];
xrmvd=[];
UPt=[];
LEt=[];
x1=1;                % for loop controlling running times of the model
% Solves the Hydrologic Balance
dt=0.01;             % (d)
tt=[0:dt:N-1];       % time axis for graphing
Nm=length(tt);
%------------------- Soil properties
soiltype=6;                  % Sand =1; Clay =11
[b,thetas, Ks1, psis,theta_w]=soil_hydraulic_values_generator(soiltype);
rhob=1500*1000;               % Soil bulk density, g/m3
Por=1.01*thetas;              % Porosity, a bit larger than soil moisture content near saturation
Ks=Ks1*1000*3600*24;          % now in mm/d
sw=0.1*(theta_w/Por);         % degree of saturation well below wilting point
s1=thetas/Por;                % soil moisture threshold for deep percolation (saturation)
s(x1,1)=0.7*s1;                 % degree of saturation at t=0, near saturation
ET=[]; LQ=[];
%------------------- Climatic condition
T=25;                      %atmospheric temperature (degree celcius)
Ca=400/1000000;             %atmospheric CO2 concentration (atm)
RH=0.54;                    %relative humidity
ET0=4;                        % Reference ET mm/d, to be modified with temperature?
VPD=0.611*exp(17.502*T/(249.91+T))*(1-RH)/101; %vapor pressure deficit, calculated with Claussius-Clapeyron equation (atm)
WUE=0.625*Ca*0.25/VPD;      %Water use efficiency at leaf level,assuming Ci/Ca=0.75 (0.37 for C4 plants), gc/gv=1/1.6 (mol CO2/mol H2O)
%------------------- Contaminant properties
kd=1e-5;                    % Equilibrium partitioning coefficient (m3/g) of contaminant adsorption
TSCF=0.5;                   % Transpiration-stream concentration factor
xo=1;
x=[];xs=[];UPx=[]; LEx=[];
x(1)=xo;                    % Contaminant concentration (g/mm3)
%------------------- Plant properties
LAI=1;                %(m2/m2)
Zr=50;               %Maximum root-zone depth (mm)
Evmax=ET0*exp(-0.398*LAI);               % surface evaporation decreases with LAI (Or & Lehmann 2019)
r=[]; UPxt=[];
k_intercept=0.05;         %Attenuation of rainwater - used for interception
UPxt(1)=0;                %Accumulated toxicant in the plant (g/mm2)
tox=0.01;                     %Parameter indicating plant sensitivity to the toxicant

for i=1:Nm
    rs=max((s(x1,i)-sw)/(s1-sw),0);   %relative saturation/effective saturation (theta_R)
    xs(i)=x(i)/(Por*s(x1,i)+rhob*kd); %solute contaminant
    dw(i)=Por*Zr;             %nZr(mm)
% Hydrologic balance
    Tr(i)=ET0*(0.33*LAI+0.45)*(-2*s(x1,i)^3+3*s(x1,i)^2)-tox*UPxt(i);            %Transpiration (mm/d);increase with LAI and s; decrease as contaminant accumulates; empirical
    Ev(i)=Evmax*rs;      %surface evaporation increases with rs (mm/d)
    ET(i)=Tr(i)+Ev(i);
    Precip=interp1(tday,Pr,tt(i),'linear');
    a_intercept=exp(-k_intercept*LAI);
    Is(i)=min(a_intercept*Precip,Ks);      %Is, ET and LQ all in mm/d
    LQ(i)=Ks*rs^(2*b+3);
    s(x1,i+1)=max(s(x1,i)+dt*(Is(i)-ET(i)-LQ(i))/dw(i),0.99*sw);
    s(x1,i+1)=min(s(x1,i+1),1);
% Contaminant mass balance    
    UPx(i)=TSCF*ET(i)*xs(i);    % Contaminant uptake (g/(mm2*d))
    LEx(i)=LQ(i)*xs(i);
    x(i+1)=max(x(i)+dt*(-UPx(i)-LEx(i))/Zr,0);  %g/mm3
    UPxt(i+1)=UPxt(i)+UPx(i)*dt;         
end
%---- Compute the efficiency at the end of period
%Eff=UPxt(Nm)/(sum(LEx*dt./Zr)+UPxt(Nm));
%effrun(nn)=UPxt(Nm+1)/(sum(LEx*dt./Zr)+UPxt(Nm+1));
%Eff(pf,st)=mean(effrun);
UPt(x1)=sum(UPx*dt);
LEt(x1)=sum(LEx*dt);
Eff(x1)=UPt(x1)/(LEt(x1)+UPt(x1));
xrmvd(x1)=1-x(Nm)/xo;
%end

%% VariablesPlottingOverTime
%----------- Plot the Hydrologic Balance Components
Plot_Soilmoisture_Time_Series
%----------- Plot the Contaminant Balance Components
Plot_Contaminant_Time_Series
%----------- Plot efficiency surface
%% plotting the track
if mod(x1,20) == 0
    figure(5)
    plot3(ones(1,Nm)*xo,s(x1,1:Nm),LAI(x1,1:Nm))
    xlabel('Initial contaminant concentration (g/mm^3)')
    ylabel('Soil moisture content')
    zlabel('LAI')
    hold on
end