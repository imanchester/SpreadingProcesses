%MAINPROBLEM2 main file to solve problem 2 with resource allocation on both
%beta and delta
% This matlab code solves Problem 2 (Risk-Constrained Resource
% Minimization) as defined in:
% https://arxiv.org/abs/2003.07555 (early version)
% https://ieeexplore.ieee.org/document/9120170 (early access)
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

% Vera Somers, June 2020, V2.0

clear all 
close all

%load csv files
Sveg=csvread('Vegetation.csv');
C=csvread('Cost2.csv');
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
Gamma=log(0.05); % risk threshold
dr = 3.3; %discount rate
betaL=1E-4; %lower bound on beta >0
deltaH=0.9; %upper bound on delta <1



%obtain state matrix 
AM=adjacencymatrix(rows,cols);
[Pveg,A]=CAmodel(Vw,theta,beta,delta,Sveg,E,AM);


%convex optimization

%parameters
y=sdpvar(n,1);
tempAM=(Pveg.*AM)';
[o,u,s] = find(tempAM);
ss=length(s);
[m,v] = size(tempAM);
rS=sdpvar(ss,1);
rij=sparse(o,u,rS,m,v);
vij=sdpvar(n,1);


Atest=A+delta*eye(n); %spreading rates only
Beta=sparse(Atest');

p0=(C/(dr*eye(n)-A'))'; %node impact vector or priority vector


betaLog=log(Beta/betaL); %upperbound on resource allocation
betaLog(isinf(betaLog)|isnan(betaLog))=0;
deltaLog=log((1-delta)/(1-deltaH));


%exponentional cone programming
Constraints= [max(log(Lambda)+y)<=Gamma,rij>=0,rij<=betaLog,vij>=0,vij<=deltaLog];

for j=1:n
    Spar=find(Beta(:,j));
    Constraints=[Constraints,logsumexp([(y(Spar)+log(Beta(Spar,j)/(1+dr))-y(j)-rij(Spar,j))' (log((1-delta)/(1+dr))-vij(j)) (log(C(j)/(1+dr))-y(j))]')<=0];
end

optimize(Constraints,sum(sum(rij))+sum(vij),sdpsettings('debug',1,'convertconvexquad',0))

KK=double(rij);
DD=double(vij);

%plot results
Plotting(Sveg,Slike,p0,KK,DD)

