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

Eff=[];
xrmvd=[];
LAIend=[];
for nn=1:18
    %nn=1;
    Eff_temp=[];
    xrmvd_temp=[];
    UPt_temp=[];
    LEt_temp=[];
    LAIend_temp=[];
    for run=1:100                % for loop controlling running times of the model
        %------------ Construct synthetic precipitation time series
        freq=1/10;                       % Return frequency between days (1/d)
        annual_precip=500;               % Annual precip. (mm/year)
        N=365;                           % Number of days to simulate the process
        dep=(annual_precip/365)/freq;
        [Pr]=precip_generate(freq,dep,N);              
        tday=[0:1:N-1];                    % Time series
        %----------- Plot synthetic precipitation time series
        %Plot_Precip_Time_Series
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
        s=[];
        s(1)=0.9*s1;                 % degree of saturation at t=0, near saturation
        ET=[]; LQ=[];
        %------------------- Climatic condition
        T=25;                      %atmospheric temperature (degree celcius)
        Ca=(400+50*(nn-1))/1000000;             %atmospheric CO2 concentration (atm)
        RH=0.5;                    %relative humidity
        ET0=4;                        % Reference ET mm/d, to be modified with temperature?
        VPD=0.611*exp(17.502*T/(249.91+T))*(1-RH)/101; %vapor pressure deficit, calculated with Claussius-Clapeyron equation (atm)
        WUE=0.625*Ca*0.25/VPD;      %Water use efficiency at leaf level,assuming Ci/Ca=0.75 (0.37 for C4 plants), gc/gv=1/1.6 (mol CO2/mol H2O)
        %------------------- Contaminant properties
        kd=1e-5;                    % Equilibrium partitioning coefficient (m3/g) of contaminant adsorption
        TSCF=0.54;                   % Transpiration-stream concentration factor
        xo=1;
        x=[];xs=[];UPx=[]; LEx=[];
        x(1)=xo;                    % Contaminant concentration (g/mm3)
        %------------------- Plant properties
        LAImax=4;                %(m2/m2)
        Zrmax=600;               %Maximum root-zone depth (mm)
        Zr=[]; r=[]; UPxt=[]; LAI=[];
        LAI(1)=1;                %Initial LAI (how much biomass has been developed before the plant is used for phytoremediation)
        %rmax=0.05;                %Maximum growth rate of LAI (1/time), noting that this value is associated with dt
        k_intercept=0.05;         %Attenuation of rainwater - used for interception
        UPxt(1)=0;                %Accumulated toxicant in the plant (g/mm2)
        tox=0.01;                     %Parameter indicating plant sensitivity to the toxicant
        SRL=120;                   %Specific root length (mm/kg)
        ma=5;                    %Empirical transfering parameter of aboveground biomass (LAI and Msh) 0.15-0.2 (kg/m2)
        Zr(1)=50;
        %------------------- Plant dynamics parameters
        kC=0.1;                    %Rate parameter of photosynthesis giving typical rate of input of substrate C (d^-1)
        kG=200;                    %Rate parameter (d^-1)
        Msh=[]; Mrt=[];    %All carbon-related variables are in per m2 area
        Msh(1)=LAI(1)*ma; Mrt(1)=Zr(1)/SRL; %Shoot biomass (kg) coupled with LAI; Root biomass (kg) coupled with Zr
        An=[]; Tr=[]; Ev=[];

        for i=1:Nm
            rs=max((s(i)-sw)/(s1-sw),0);   %relative saturation/effective saturation (theta_R)
            xs(i)=x(i)/(Por*s(i)+rhob*kd); %solute contaminant
            dw(i)=Por*Zr(i);             %nZr(mm)
        % Plant dynamics
            Tr(i)=ET0*(0.33*LAI(i)+0.45)*(-2*s(i)^3+3*s(i)^2)-tox*UPxt(i);            %Transpiration (mm/d);increase with LAI and s; decrease as contaminant accumulates; empirical
            An(i)=12*Tr(i)*WUE/18;               %Photosynthesis rate kg C/(m2*d); note that there is a unit conversion (mm/d to kg/(m2*d)); 12 and 18 are MW of C and H20       
        % Hydrologic balance
            Evmax=(ET0-Tr(i))*exp(-0.398*LAI(i));               % surface evaporation decreases with LAI (Or & Lehmann 2019)
            Ev(i)=Evmax*rs;      %surface evaporation increases with rs (mm/d)
            ET(i)=Tr(i)+Ev(i);
            Precip=interp1(tday,Pr,tt(i),'linear');
            a_intercept=exp(-k_intercept*LAI(i));
            Is(i)=min(a_intercept*Precip,Ks);      %Is, ET and LQ all in mm/d
            LQ(i)=Ks*rs^(2*b+3);
            s(i+1)=max(s(i)+dt*(Is(i)-ET(i)-LQ(i))/dw(i),0.99*sw);
            s(i+1)=min(s(i+1),1);
        % Contaminant mass balance    
            UPx(i)=TSCF*ET(i)*xs(i);    % Contaminant uptake (g/(mm2*d))
            LEx(i)=LQ(i)*xs(i);
            x(i+1)=max(x(i)+dt*(-UPx(i)-LEx(i))/Zr(i),0);  %g/mm3
            UPxt(i+1)=UPxt(i)+UPx(i)*dt;         
        % Carbon assimilation and partitioning
             Time_2_max_LAI=3650; %days needed to reach the maximum LAI
             Re(i)=(Msh(i)+Mrt(i))/Time_2_max_LAI;%An(i)/2;                 %Respiration counts for half of total C assimilation. Farrar 1985, Amthor 1989
             dMsh=max(0.75*(An(i)-Re(i))*2*dt,0);         %Plant tissue typically contains about 45-50% carbon (so (An-Re)*2); Assuming the biomass partition (0.75 aboveground) follows the empirical beta from Niklas 2005
             Msh(i+1)=Msh(i)+dMsh;
             LAI(i+1)=min(Msh(i+1)/ma,LAImax);
             dMrt=max(0.25*(An(i)-Re(i))*2*dt,0);
             Mrt(i+1)=Mrt(i)+dMrt;
             Zr(i+1)=min(Mrt(i+1)*SRL,Zrmax);
        end
        UPt_temp(run)=sum(UPx*dt);
        LEt_temp(run)=sum(LEx*dt);
        Eff_temp(run)=UPt_temp(run)/(LEt_temp(run)+UPt_temp(run));
        xrmvd_temp(run)=1-x(Nm)/xo;
        LAIend_temp(run)=LAI(Nm);
    end
    %---- Compute the efficiency at the end of period
    %Eff=UPxt(Nm)/(sum(LEx*dt./Zr)+UPxt(Nm));
    %effrun(nn)=UPxt(Nm+1)/(sum(LEx*dt./Zr)+UPxt(Nm+1));
    %Eff(pf,st)=mean(effrun);
    Eff(nn)=mean(Eff_temp);
    xrmvd(nn)=mean(xrmvd_temp);
    LAIend(nn)=mean(LAIend_temp);
end

%% Variables plotting over time for scenarios
figure(6)
subplot(3,1,1)
plot((400:50:1250),Eff,'k-')
ylabel ('Efficiency','fontweight','normal','fontsize',15)

subplot(3,1,2)
plot((400:50:1250),xrmvd,'k-')
%xlabel('Annual rainfall (mm)','fontweight','normal','fontsize',15)
%xlabel('Rainfall return period (d)','fontweight','normal','fontsize',15)
ylabel ('Removed x %','fontweight','normal','fontsize',15)

subplot(3,1,3)
plot((400:50:1250),LAIend,'k-')
ylabel('Final LAI','fontweight','normal','fontsize',15)
xlabel('CO_2 concentration (ppm)','fontweight','normal','fontsize',15)
%xlabel ('Temperature (\circC)','fontweight','normal','fontsize',15)
%% VariablesPlottingOverTime
%----------- Plot the Hydrologic Balance Components
Plot_Soilmoisture_Time_Series
%----------- Plot the Contaminant Balance Components
Plot_Contaminant_Time_Series
%----------- Plot the plant dynamics
Plot_plant_dynamics
%----------- Plot efficiency surface
%% plotting the track
if mod(nn,100) == 0
    figure(5)
    plot3(ones(1,Nm)*xo,s(x1,1:Nm),LAI(x1,1:Nm))
    xlabel('Initial contaminant concentration (g/mm^3)')
    ylabel('Soil moisture content')
    zlabel('LAI')
    hold on
end