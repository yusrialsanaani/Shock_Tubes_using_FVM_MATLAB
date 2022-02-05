%% =========================================================================
%  Matlab program to solve Shock Tubes problem.
%  Done by Yusri Al-Sanaani
%% Initialization
clear all;  close all;    clc
L=1.5;      N=80;         dx=L/(N-1);     x=0:dx:L;
tmax=0.3;   dt=0.0001;    sigma=-dt/dx;   gamma=1.4;  t=0;
%% Define the initial left & right parameters of Sod's shock tube,[rho,u,p]
rhoR(1:N/2)=0.125;     uR(1:N/2)=0;      pR(1:N/2)=0.1;
rhoL(1:N/2)=1;         uL(1:N/2)=0;      pL(1:N/2)=1;
%% Define the left & right variables of W=[rho, rho*u, rho*E] based on ui=[rho,u, p]
wR(1,:)=rhoR; wR(2,:)=rhoR.*uR; rhoER=(pR/(gamma-1))+0.5*wR(2,:).*uR; wR(3,:)=rhoER;    
wL(1,:)=rhoL; wL(2,:)=rhoL.*uL; rhoEL=(pL/(gamma-1))+0.5*wL(2,:).*uL; wL(3,:)=rhoEL;
w=[wL wR];
%% Solve the first-order upwind Roe’s scheme using Forward Euler Methods
while t<=tmax
    t=t+dt;
    F=RoeFlux(w,gamma,N);
   for i=2:N-1
       w(:,i)=w(:,i)+sigma*(F(:,i)-F(:,i-1));
   end
end
%% obtain the rho, u, p, and energy from final W
rho=w(1,:); u=w(2,:)./rho; p=(gamma-1)*(w(3,:)-0.5*w(2,:).*u); E=w(3,:)./rho;
%--------------------------------------------------------------------------
figure
plot(x,rho,'LineWidth',2.0);grid off 
title({['Density at \itt = ' num2str(t),' and N = ',num2str(N)]})
xlabel('x','fontweight','bold');ylabel('\rho','fontweight','bold');
set(gca,'FontName', 'Times New Roman','FontSize',12,'linewidth',1.5,'box','off','XAxisLocation','origin');
%--------------------------------------------------------------------------
figure
plot(x,u,'LineWidth',2.0);
title({['Velocity at \itt = ' num2str(t),' and N = ',num2str(N)]})
xlabel('x','fontweight','bold');ylabel('u','fontweight','bold');
set(gca,'FontName', 'Times New Roman','FontSize',12,'linewidth',1.5,'box','off','XAxisLocation','origin');
%--------------------------------------------------------------------------
figure
plot(x,p,'LineWidth',2.0);
title({['Pressure at \itt = ' num2str(t),' and N = ',num2str(N)]})
xlabel('x','fontweight','bold');ylabel('p','fontweight','bold');
set(gca,'FontName', 'Times New Roman','FontSize',12,'linewidth',1.5,'box','off','XAxisLocation','origin');
%--------------------------------------------------------------------------
figure
plot(x,E,'LineWidth',2.0);
title({['Energy at \itt = ' num2str(t),' and N = ',num2str(N)]})
xlabel('x','fontweight','bold');ylabel('E','fontweight','bold');
set(gca,'FontName', 'Times New Roman','FontSize',12,'linewidth',1.5,'box','off','XAxisLocation','origin');
%%


