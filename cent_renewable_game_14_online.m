clc;clear all;
% #Players = 14
tic

%% Load data
load YBUS_14.csv
load Gen_coeffs.csv

%% Initialization
n=14;                                                                       %Total agents
p=3;                                                                        %Slack bus no.
T=1800;                                                                     %Horizon

%Case-1
p_g = zeros(n,1);                                                           %Generation(MW)      
p_l = zeros(n,1);                                                           %P-Load (MW)
p_l(3,1)=95; p_l(4,1)=125; p_l(6,1)=100; p_l(8,1)=106; p_l(10,1)=115;     
p_l(12,1)=90; p_l(14,1)=80;
p_gmax=[280;110;0;0;0;0;70;160;0;0;0;130;0;0];                              %Generation capacity (MW)  
p_g = p_g/100;                                                              %PU base = 100MVA
p_l = p_l/100;
p_gmax = p_gmax/100;

%Case-2
a = [0.001 + (0.08-0.001).*rand(5,T)];                     %Thermal generation cost parameters  
b = [1 + (5-1).*rand(5,T)];
c = [zeros(5,T)];

%% Topology (IEEE-14 bus)
Y = YBUS_14;

%% Centralized Optimization
for k=1:1
    cvx_begin
    variables pg(n) theta(n)   
    expressions cost
    j=1;
    for i=1:n
        if i==1 || i==2 || i==7 || i==8 || i==12
            cost = cost + (a(j,k)*pg(i)^2 + b(j,k)*pg(i) + c(j,k));
            j=j+1;
        end
    end
    minimize cost
    subject to 
    theta(p)==0;                                                          %Slack Bus
    for i=1:n
        if i~=p
            pg(i)==p_l(i)+theta(i)*sum(Y(i,:))-Y(i,:)*theta;
            pg(i)>=0;
            pg(i)<=p_gmax(i);
        end
    end
    pg(p)==p_l(p)+theta(p)*sum(Y(p,:))-Y(p,:)*theta;
    pg(p)>=0;
    pg(p)<=p_gmax(p);
    cvx_end 
end
cent_pg = 100*pg;
cent_theta = theta;
















