%%%%%
%Online Optimal Power Flow
%%%%%

cent_renewable_game_14_online
clc;
clearvars -except p_gmax p_l a b c n p T Y
% #Players = 14
tic

%% Initialization
mu = ones(n,1);                                                                 %Dual variable
mu_tilde = zeros(n,1);                                                          %Weighted Dual variable
theta = zeros(n,1);                                                             %Load angle                                                               
cons_reg = [0];                                                                 %Avg. Equality Constraint Violation Regret

p_g = zeros(n,1);                                                               %Generation      
p_g = p_g/100;                                                                  %PU base = 100MVA
reg = zeros(n,T);
avg_reg = zeros(n,T);
costo = zeros(n,1);
viol = zeros(n,1);
violcum = zeros(n,T);
avg_violcum = zeros(n,T);

%% Weighted Adjacency Matrix Generation
TP = Y~=0;
cvx_begin quiet
variable W(n,n) 
minimize 0
subject to
for i=1:n
    for j=1:n
        if TP(i,j)~=1 
            W(i,j)==0;
        end
        W(i,j)>=0;
    end
end
W*ones(n,1)==ones(n,1);
ones(1,n)*W==ones(1,n);
cvx_end

%% Distributed Optimization
for k=1:T 
    %Dual Update (mu_tilde)
    mu_tilde(:,k) = W*mu(:,k);
    %Primal Update (p_g)
    kj=1;
    for i=1:n
        if i==1 || i==2 || i==7 || i==8 || i==12                                %Slack bus and Generators
            P = p_g(i,k) - (1/sqrt(k))*(2*a(kj,k)*p_g(i,k) + b(kj,k) + mu_tilde(i,k));
            kj=kj+1;
            if P >= p_gmax(i)
                p_g(i,k+1) = p_gmax(i);
            elseif P <= 0
                p_g(i,k+1) = 0;
            else
                p_g(i,k+1) = P;
            end
        else                                                                    %Loads 
            p_g(i,k+1)=p_g(i,k);
        end
    end
    %Primal Update (theta)
    for i=1:n
        th=0;
        if i==p                                                                 %Slack bus
            theta(i,k+1) = 0;
        else                                                                    %Thermal Generators and Loads     
            th = theta(i,k) - (1/sqrt(k))*(p_g(i,k)-p_l(i)-theta(i,k)*sum(Y(:,i))+Y(i,:)*theta(:,k));
            if th>=pi/2
                theta(i,k+1) = pi/2;
            elseif th<=-pi/2
                theta(i,k+1) = -pi/2;
            else
                theta(i,k+1) = th;
            end
        end
    end
    %Dual Update (mu)
    for i=1:n
        mu(i,k+1) = mu_tilde(i,k) + (1/sqrt(k))*(p_g(i,k)-p_l(i)-theta(i,k)*sum(Y(:,i))+Y(i,:)*theta(:,k))... 
        - (1/k)*(mu_tilde(:,k)'*Y(:,i) - mu_tilde(i,k)*sum(Y(:,i)));
    end
    
    %% Static Regret
    %Optimal Decision in hindsight
    cvx_begin quiet
    variables pg(n) thetap(n)   
    expressions cost
    for ki=1:k
        ji=1;
        for i=1:n
            if i==1 || i==2 || i==7 || i==8 || i==12
                cost = cost + (a(ji,ki)*pg(i)^2 + b(ji,ki)*pg(i) + c(ji,ki));
                ji=ji+1;
            end
        end
    end
    minimize cost
    subject to 
    thetap(p)==0;                                                               %Slack Bus
    for i=1:n
        if i~=p
            pg(i)==p_l(i)+thetap(i)*sum(Y(i,:))-Y(i,:)*thetap;
            pg(i)>=0;
            pg(i)<=p_gmax(i);
        end
    end
    pg(p)==p_l(p)+thetap(p)*sum(Y(p,:))-Y(p,:)*thetap;
    pg(p)>=0;
    pg(p)<=p_gmax(p);
    cvx_end
    jki=1;
    costc = zeros(n,1);
    for i=1:n
        if i==1 || i==2 || i==7 || i==8 || i==12
            for kj=1:k
                costc(i) = costc(i) + a(jki,kj)*pg(i)^2 + b(jki,kj)*pg(i) + c(jki,kj);
            end
            jki=jki+1;
        end
    end
    %Online Decision
    jo=1;
    for i=1:n
        if i==1 || i==2 || i==7 || i==8 || i==12
            costo(i) = costo(i) + (a(jo,k)*p_g(i,k)^2 + b(jo,k)*p_g(i,k) + c(jo,k));
            jo=jo+1;
        end
    end
    %Individual Regret calculation
    for i=1:n
        if i==1 || i==2 || i==7 || i==8 || i==12
            reg(i,k) = costo(i) - costc(i);
            avg_reg(i,k) = (1/k)*reg(i,k);
        end
    end
    
    %% Equality Constraint Violation 
    for i=1:n
        viol(i) = viol(i) + (p_g(i,k)-p_l(i)-theta(i,k)*sum(Y(:,i))+Y(i,:)*theta(:,k));
        violcum(i,k) = viol(i);
        avg_violcum(i,k) = (1/k)*violcum(i,k);
    end
    toc   
    [k, sum(avg_reg(:,k)), sum(avg_violcum(:,k)), sum(reg(:,k))/100, sum(violcum(:,k))/100]
end

%% Performance plot
figure(1)
plot([1:k],sum(avg_reg,1),[1:k],sum(avg_violcum,1),'LineWidth',0.8);               %N/W Performance (Averaged over time)                                       
title('Performance (averaged over T)');
xlabel(['T'],'interpreter','latex','FontWeight','bold','FontSize',15);
ylabel(['$\frac{R(T)}{T}$'],'interpreter','latex','FontWeight','bold','FontSize',20);
legend('N/W Regret','N/W Violation');

figure(2)
plot([1:k],sum(reg,1),[1:k],sum(violcum,1),'LineWidth',0.8);                       %N/W Performance                                       
title('Performance');
xlabel(['T'],'interpreter','latex','FontWeight','bold','FontSize',15);
ylabel(['$R(T)$'],'interpreter','latex','FontWeight','bold','FontSize',20);
legend('N/W Regret','N/W Violation');

figure(3)
subplot(2,1,1);
for i=1:n
    if i==1 || i==2 || i==7 || i==8 || i==12
        plot([1:k],avg_reg(i,:),'LineWidth',0.8);                                  %Agent Regret (Averaged over time)   
        hold on
    end
end
title('Agents Regret (averaged over T)');
xlabel(['T'],'interpreter','latex','FontWeight','bold','FontSize',15);
ylabel(['$\frac{R_i(T)}{T}$'],'interpreter','latex','FontWeight','bold','FontSize',20);
subplot(2,1,2);
for i=1:n
    if i==1 || i==2 || i==7 || i==8 || i==12
        plot([1:k],avg_violcum(i,:),'LineWidth',0.8);                              %Agent Violation (Averaged over time)  
        hold on
    end
end
title('Agents Violation (averaged over T)');
xlabel(['T'],'interpreter','latex','FontWeight','bold','FontSize',15);
ylabel(['$\frac{R^{ec}_i(T)}{T}$'],'interpreter','latex','FontWeight','bold','FontSize',20);

figure(4)
subplot(2,1,1);
for i=1:n
    if i==1 || i==2 || i==7 || i==8 || i==12
        plot([1:k],reg(i,:),'LineWidth',0.8);                                      %Agent Regret   
        hold on
    end
end
title('Agents Regret');
xlabel(['T'],'interpreter','latex','FontWeight','bold','FontSize',15);
ylabel(['$R_i(T)$'],'interpreter','latex','FontWeight','bold','FontSize',20);
subplot(2,1,2);
for i=1:n
    if i==1 || i==2 || i==7 || i==8 || i==12
        plot([1:k],violcum(i,:),'LineWidth',0.8);                                  %Agent Violation   
        hold on
    end
end
title('Agents Violation');
xlabel(['T'],'interpreter','latex','FontWeight','bold','FontSize',15);
ylabel(['$R^{ec}_i(T)$'],'interpreter','latex','FontWeight','bold','FontSize',20);

figure(5)
for i=1:n
    p5=plot([1:k],mu(i,2:end),'.--','LineWidth',0.4);                              %LMP                         
    hold on
end
title('Locational Marginal Prices');
xlabel(['k'],'interpreter','latex','FontWeight','bold','FontSize',15);
ylabel(['$\mu(k)$'],'interpreter','latex','FontWeight','bold','FontSize',15);

figure(6)
for i=1:n
    if i==1 || i==2 || i==7 || i==8 || i==12 
        p1=plot([1:k],100*p_g(i,2:end),'.--','LineWidth',1.2);                     %Distributed Power
        hold on
    end
end
title('Optimal Generated Power (p_g)');
xlabel(['k'],'interpreter','latex','FontWeight','bold','FontSize',20);
ylabel(['$p_g(k)$'],'interpreter','latex','FontWeight','bold','FontSize',20);
legend('DM 1','DM 2','DM 7','DM 8','DM 12');

figure(7)
for i=1:n
    if i==1 || i==2 || i==7 || i==8 || i==12 || i==p
        p3=plot([1:k],theta(i,2:end),'.--','LineWidth',0.1);                       %Distributed Theta 
        hold on
    end
end
title('Optimal Load Angle (\theta)');
xlabel(['k'],'interpreter','latex','FontWeight','bold','FontSize',20);
ylabel(['$\theta(k)$'],'interpreter','latex','FontWeight','bold','FontSize',20);
legend('DM 1','DM 2','DM 7','DM 8','DM 12');
