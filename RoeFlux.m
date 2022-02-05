function F=RoeFlux(w,gamma,N)
rho=w(1,:);                            % Get the densities                
u=w(2,:)./w(1,:);                      % Get the velocities
p=(gamma-1)*(w(3,:)-0.5*w(2,:).*u);    % Get the pressures
E=w(3,:)./rho;                         % Get the specific energy
H=E+p./rho;                            % Get the Enthalpy, H=E+p/rho
% Get the flux
F_w(1,:)= w(2,:); 
F_w(2,:)= w(2,:).*u+p;
F_w(3,:)= (w(3,:)+p).*u;
%% Calculate Roe average values
for j=1:N-1             
  u_avg=(u(:,j)*sqrt(rho(:,j))+u(:,j+1)*sqrt(rho(:,j+1)))/(sqrt(rho(:,j))+sqrt(rho(:,j+1))); % Roe's Velocity
  H_avg=(H(j)*sqrt(rho(:,j))+H(j+1)*sqrt(rho(:,j+1)))/(sqrt(rho(:,j))+sqrt(rho(:,j+1)));     % Roe's Enthalpy
  a=sqrt((gamma-1.0)*(H_avg-0.5*u_avg^2));  % the speed of sound
  beta=1/(2*a^2);
  R=[1 beta beta;u_avg (u_avg+a)*beta (u_avg-a)*beta;0.5*u_avg^2 beta*(H_avg+a*u_avg) beta*(H_avg-a*u_avg)];
  L=inv(R);
  D=diag([abs(u_avg) abs(u_avg+a) abs(u_avg-a)]);
  %Apply Haartens entropy fix
del=0.001;
  for i=1:3
      if(i==1||i==3)
          if abs(D(i,i))<del
              D(i,i)=0.5*(D(i,i)^2/del+del);
          end
      end
  end
  wL=w(:,j+1);wR=w(:,j);
  % Calculate F=0.5*(F(WL}+F(WR))-R*D*L*(WL-WR)
  F(:,j)=0.5*(F_w(:,j)+F_w(:,j+1))-0.5*R*D*L*(wL-wR);
end
end



