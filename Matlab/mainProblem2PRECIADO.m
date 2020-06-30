%MAINPROBLEM2 main file to solve problem 2 with resource allocation on beta only 
% This matlab code solves Problem 2 (Resource-Constrained Risk
% Minimization) using our proposed resource model while minimizing the
% dominant eigenvalue (as done in Preciado et al, 2014)
% 
% USE MOSEK 8 (or older versions to be able to solve with GP programming)
%
% Inputs:
% - csv files containing vegetation, cost, likelihood and elevation
% data
% - wind speed and direction
% - baseline spreading rate and recovery rate
% - lower bound on beta
% - discount rate
% - risk threshold
%
% Outputs:
% - resource allocation map
% - risk map

% Vera Somers, March 2020

clear all 
close all

%load csv files
Sveg=csvread('Vegetation.csv');
C=csvread('Cost.csv');
C=C(:)';
Slike=csvread('Likelihood.csv');
E=csvread('Elevation.csv');


%parameter set-up
[rows,cols]=size(Slike); %define number of rows and columns in grid based graph
n=rows*cols; %define number of nodes in graph
Lambda=Slike(:); %likelihood

%decision variables
Vw=4; %wind speed (m/s)
theta= 225; %wind direction (degrees)
delta=0.2; %recovery rate
beta=0.5; %infection rate base line
Gamma=0.0516; % risk threshold
dr = 3.5; %discount rate
betaL=1E-4; %lower bound on beta


%obtain state matrix 
AM=adjacencymatrix(rows,cols);
[Pveg,A]=CAmodel(Vw,theta,beta,delta,Sveg,E,AM);


%convex optimization

%parameters
p=sdpvar(n,1);
tempAM=(Pveg.*AM)';
[o,u,s] = find(tempAM);
ss=length(s);
[m,v] = size(tempAM);
betaS=sdpvar(ss,1);
BetaN=sparse(o,u,betaS,m,v);


Atest=A+delta*eye(n); %spreading rates only
Beta=Atest';

p0=(C/(dr*eye(n)-A'))'; %node impact vector or priority vector


betaH=nonzeros(Beta(:)); %upperbound on spreading rate
betaMax=max(Beta(:));

%exponentional cone programming
Constraints= [p>=0,betaS>=betaL,BetaN<=Beta,max(Lambda.*p)<=Gamma];

for j=1:n
    Constraints=[Constraints,sum(p.*BetaN(:,j))/(p(j)*delta) - dr/delta + C(j)/(p(j)*delta) <=1];
end

optimize(Constraints,sum((1-betaS(:)./betaH(:))./(betaS(:)/betaL - betaS(:)/betaMax)),sdpsettings('solver','mosek-geometric','debug',1,'convertconvexquad',0))

KK=-log(double(BetaN)./Beta);

KK(isinf(KK)|isnan(KK))=0;

%plot results
Plotting(Sveg,Slike,p0,KK)













