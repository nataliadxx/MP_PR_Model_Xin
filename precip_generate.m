%--------------------------------------------------------------------------
%precip_generate.m generates randomly daily rainfall with 
%depth and interpulse durations from exponential distribution. 
%The daily precipitation approximates a marked Poisson process
%https://en.wikipedia.org/wiki/Point_process
%
%Pr = daily rainfall
%dep = mean rainfall depth (mm/d) 
%freq = return frequency of rainfall (1/d) 
%N = length of series to be generated (days)
%--------------------------------------------------------------------------
function [Pr]=precip_generate (freq,dep, N)
P=[];Pr=[];
Pr=zeros(1,N);
t=[1:1:N];
%--------- Generate Durations and Depth of Precipitation
r=rand(1,N);    %vector of random numbers [0,1] N-long.
%----- Generate exponentially distributed period between rain events
beta=freq;
tau=(-1/beta)*(log(1-r+eps));   
%----- Generate exponentially distributed rainfall depth
beta=1/dep;
eta=(-1/beta)*(log(1-r+eps));   
%----- Construct vector of times when it rained
tt=cumsum(tau)-tau(1)+1;
tti=floor(tt);
P(tti)=eta;
Pr=P(1:N); %vector of rain depths of length N


