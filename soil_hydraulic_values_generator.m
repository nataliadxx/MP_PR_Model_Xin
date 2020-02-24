function [b,thetas, Ks, psis,thetaw]=soil_hydraulic_values_generator(m)
%	---------------------- Variables/Units-------------------------------------
%	b       = Empirical Coefficient that varies with soil texture
%	Ks      = Reference hydraulic conductivity close to saturation (m/s)
%	Psis    = Reference soil pressure close to saturation (m)
%   thetas  = Soil moisture content close to saturation
%	theta_w = Wilting point
%	Soils 1 [sand] to 11 [clay] are from Clapp and Hornberger (1978).  
%	Soil 12 is measured at Duke Forest
%   m       = soil number (1-->12)
%   Parameters used in soil water characteristic curves and hydraulic conductivity
%   Psi=Psis(theta/thetas)^(-b)
%   K = Ks (theta/thetas)^(2b+3)
%-----------------------------------------------------------------------------
	b1(1)=4.05; %sand
	b1(2)=4.38; %loamy sand
	b1(3)=4.90; %sandy loam
	b1(4)=5.30; %silt loam
	b1(5)=5.39; %loam
	b1(6)=7.12; %sandy clay loam
	b1(7)=7.75; %silty clay loam
	b1(8)=8.52; %clay loam
	b1(9)=10.4; %sandy clay
	b1(10)=10.4; %silty clay
	b1(11)=11.4; %clay
	b1(12)=2.79;
	
	thetas1(1)=0.395; 
	thetas1(2)=0.441; 
	thetas1(3)=0.435; 
	thetas1(4)=0.485; 
	thetas1(5)=0.451; 
	thetas1(6)=0.42; 
	thetas1(7)=0.477; 
	thetas1(8)=0.476; 
	thetas1(9)=0.426; 
	thetas1(10)=0.492; 
	thetas1(11)=0.482;
	thetas1(12)=0.54;
   

	Ks1(1)=176e-6;
	Ks1(2)=156e-6;
	Ks1(3)=34.1e-6; 
	Ks1(4)=7.2e-6 ;
	Ks1(5)=7.0e-6 ;
	Ks1(6)=6.3e-6 ;
	Ks1(7)=1.70e-6 ;
	Ks1(8)=2.5e-6 ;
	Ks1(9)=2.2e-6 ;
	Ks1(10)=1.0e-6 ;
	Ks1(11)=1.3e-6;
	Ks1(12)=0.9269e-6;
   	                               
	Psis1(1)=-0.121; 
	Psis1(2)=-0.09; 
	Psis1(3)=-0.218; 
	Psis1(4)=-0.786; 
	Psis1(5)=-0.478; 
	Psis1(6)=-0.299; 
	Psis1(7)=-0.356; 
	Psis1(8)=-0.63; 
	Psis1(9)=-0.153; 
	Psis1(10)=-0.490; 
	Psis1(11)=-0.405;                
	Psis1(12)=-0.080;
   
	
	theta_w1(1)=0.07; 
	theta_w1(2)=0.075; 
	theta_w1(3)=0.114; 
	theta_w1(4)=0.1794; 
	theta_w1(5)=0.1547; 
	theta_w1(6)=0.1749; 
	theta_w1(7)=0.2181; 
	theta_w1(8)=0.2498; 
	theta_w1(9)=0.2193; 
	theta_w1(10)=0.2838; 
	theta_w1(11)=0.2864;
	theta_w1(12)=0.10;
   
	Ks=Ks1(m);
	b=b1(m);
	psis=Psis1(m);
	thetas=thetas1(m);
    thetaw=theta_w1(m);
