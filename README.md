# multi-tissue-viscous-model

## Run the simulation:
To run the simulation you can run the following command in your Matlab command window:

multi_tissue_brinkman(N,filename,Tmax,g1,g2,injpsm,injnt,beta1,beta2,mu,e,thresh_pz,alpha,nmax_growth, nmax)

With the parameter names replaced with their values (for example the parameters listed below for the wild-type case)

## Parameters for WT simulations 

The values of the parameters that can be tested to obtain the wild-type simulations :

nmax_growth=0.9*nmax;  %maximal density for growth before stopping growth
Tmax=272000 %20h max simulation time in seconds
beta1=0.2*10^5; %PSM viscosity
beta2=0.3*10^5; %NT viscosity
e=3; %pressure sensitivity
mu=10^(12); %friction
g1=1/31500; %s^-1 psm growth rate
g2=1/38988; %s^-1 nt growth rate
injpsm= 6.405040788/3600  %psm injection 6.4 cell/h 
injnt=4.67559029/3600 % 4.6 nt injection cell/h
nmax= 0.02/10^(-12);    %maximal density reached
alpha = 8; %steepness of the growth function 
N=21; %number of mesh cells 
thresh_pz=10^(10); %the density to reach for the PZ tip to evolve
filename='choose_your_favorite_filename';

see supplementary materials and methods for more details

## PSM and notochord
Use the same code to simulate the PSM and Notochord dynamics; only change the parameters for the notochord (see supplementary materials and methods)


## Obtain the simulation video

To obtain the simulation video run the code generate_video_simulation.m 
In this code all you need to do is replace the filename variable with the filename you gave when  you ran multi_tissue_brinkman.m

## Obtain simulation analysis and plots

Simply run the code simulation_analysis.m
