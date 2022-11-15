clc
clear all
close all

dt = 0.002; % time step
dy = 0.001; % Grid spacing
ly = 0.04;  % Distance between the walls
n = ly/dy;  % Number of intervals

nstep = 540; % Number of time step


U0 = 40; % Velocity of lower wall

nu = 0.000217; % Kinematic visocsity of oil

% grid generation 

   for i = 1:n+1
    y(i) = (i-1)*dy;
   end




% Initialization(velocity at t = 0)
for j = 2:n 
    u(1,j) = 0;
end


% Boundary conditions
for i = 1:nstep+1
    for j = 1:n+1
        if j == 1
            u(i,j) = U0;
        end
        if j == n+1
            u(i,j) = 0;
        end
    end
end

p = (2*nu*dt)/(dy)^2;

for j = 1:n-1
 a(j) = p/2;
 b(j) = -(p+1);
 c(j)= a(j);
end 

% Applying Laasonen implicit to find initial set of data
% Create tridiagonal matrix 

 P = diag(a(1:n-2),-1)+diag(b(1:n-1))+diag(c(1:n-2),1);
 
 
 for i = 1:n-1
     if i==1
     Q(i) = -u(1,i+1)-a(i)*u(2,1);
     elseif i==n-1
     Q(i) = -u(1,i+1)-c(i)*u(2,n+1); 
     else
      Q(i) = -u(1,i+1);
     end
 end
 
 R = P\Q';
 
 
% Code part of Discretized equation 3.11           
for i = 1:nstep 
  for j = 2:n
      if i==1
       u(i+1,j) = R(j-1);
      else 
        u(i+1,j) = ((1-p)*u(i-1,j)+p*(u(i,j+1)+u(i,j-1)))/(1+p);
      end
  end
end

% creating file and writing data on it
 fid = fopen( 'couette_flow_DuFort1.txt', 'wt' );
  fprintf(fid,'       y        t = 0    t = 0.18   t = 0.36   t = 0.54   t = 0.72   t = 0.90   t = 1.08 \n');
  for i = 1:n+1
  A = [y(i);u(1,i);u(91,i);u(181,i);u(271,i);u(361,i);u(451,i);u(541,i)];
   fprintf(fid,'%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n',A);
  end
  
fclose(fid)

% plotting velocity profile at different time steps
figure(1)
plot(u(1,:),y,'-*r',u(91,:),y,'-or',u(181,:),y,'-sr',u(271,:),y,'-xr',u(361,:),y,'-+r',u(451,:),y,'-pr',u(541,:),y,'-dr');
legend('T=0.00 sec','T=0.18 sec','T=0.36 sec','T=0.54 sec','T=0.72 sec','T=0.90 sec','T=1.08 sec')
 xlabel('u (m/s)','fontsize', 15)
 ylabel('y (m)','fontsize', 15)


% Plotting Syntax to animate the profile 
% at different time steps
figure(2)    
for i = 1:nstep
 
        h = plot(u(i,:),y(:),'.r');
        set(gca, 'fontsize',15)
        axis([0 42 0 ly])
        xlabel('u','fontsize', 15)
        ylabel('y','fontsize', 15)
        set(h, 'XData',u(i,:),'YData',y(:));
        drawnow 
   
end

