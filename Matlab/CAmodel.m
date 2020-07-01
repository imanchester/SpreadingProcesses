function [Pveg,A]=CAmodel(Vw,theta,beta,delta,Sveg,E,AM)
%CAMODEL sparse state matrix for spreading processes 
% This function creates the state matrix A consisting of the spreading and
% recovery rates for a spreading process based on vegetation, wind and
% elevation. Total spreading rate
% beta=beta_base*beta_veg*beta_wind*beta_elevation, based on the CA models
% presented by:
% I. Karafyllidis and A. Thanailakis, “A model for predicting forest fire spreading using cellular automata,” 1997
% A. Alexandridis, D. Vakalis, C. I. Siettos, and G. V. Bafas, “A cellular
% automata model for forest fire spread prediction: The case of the wildfire that swept through Spetses Island in 1990,” 2008.
%
% Inputs:
% - Vw (wind speed (m/s))
% - theta (wind direction (degrees))
% - beta (baseline spreading rate)
% - delta (baseline recovery rate)
% - Sveg (vegetation state matrix, eucalyptic forest=1, grassland=2,
% Desert=3, city=4, water/unburnable=5)
% - E (elevation (m))
% - AM (adjacency matrix)
%
% Outputs:
% - Pveg (spreading rate correction based on vegetation)
% - A (state matrix)

% Vera Somers, June 2020, V2.0 (updated city spreading)

%parameter set-up, following cited papers
c1=0.045;
c2=0.131;
s_a=0.078; %slope correction factor
g=10; %grid size
[rows,cols]=size(Sveg);
n=rows*cols;

%setting up beta_veg correction factors
Pveg=ones(rows,cols);

for i1=1:rows
    for j1=1:cols
        if Sveg(i1,j1)==1
           Pveg(i1,j1)=1.4;
        elseif Sveg(i1,j1)==3
           Pveg(i1,j1)=0.1;
        elseif Sveg(i1,j1)==4
            Pveg(i1,j1)=0.5;  
        elseif Sveg(i1,j1)==5
            Pveg(i1,j1)=0;
        end
    end
end

Pvegtemp=Pveg; %to correct resource allocation city in mainProblem (new in V2.0)

for i1=1:rows
    for j1=1:cols
        if Sveg(i1,j1)==4
           Pvegtemp(i1,j1)=0;  
        end
    end
end

Pveg=Pveg(:);   
Pveg=repmat(Pveg,1,n);
Pvegtemp=Pvegtemp(:);
Pvegtemp=repmat(Pvegtemp,1,n);


E=E(:);

rX=repmat(linspace(0,cols,cols),[rows 1]);
rX=rX(:)';
rY=repmat(linspace(rows,0,rows),[1 cols]);
PW=zeros(n);
PS=zeros(n);

%creating PW=beta_ij based on wind, vegetation and elevation
for i=1:n
    for j=1:n
            theta_ij=atan2d(rY(j)-rY(i),rX(j)-rX(i));
            theta_w=abs(theta_ij-theta);
            pveg=Pveg(i,j);
            pw = beta*pveg*exp(c1*Vw)*exp(Vw*c2*(cosd(theta_w)-1));
            theta_s=atan((E(i)-E(j))/g); %adjacent
            if rX(i)~=rX(j) && rY(i)~=rY(j)
               theta_s=atan((E(i)-E(j))/(g*sqrt(2))); %diagonal
               pw=pw*2*(sqrt(2)-1);
            end
            ps=exp(s_a*theta_s);
            pw=ps*pw;
            if pw > 1
                pw = 1;
            elseif pw < 0
                pw =0; 
            end
            PW(i,j)=pw;
            PS(i,j)=ps;
    end
end

%final state matrix
if Vw==0 %no wind
    A=(-delta*eye(n)+beta.*Pveg.*PS.*AM);   
else
    A=(-delta*eye(n)+PW.*AM);
end


Pveg=Pvegtemp;


end


