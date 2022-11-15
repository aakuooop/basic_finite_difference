clc
clear all
close all

% Given input data


a = 300;                                          % speed of sound at the altitude of 10 km
dx =5.0;                                          % Grid size
dt = 0.01666;                                     % step size
c = a*dt/dx;                                      % Courant number
lx = 300;                                         % length of one-dimensional tube
m = lx/dx;                                        % No. of intervals
nstep = 120;                                      % No. of time steps
U1 = 0;
U2 = 0;
% Grid generation

for i = 1:m+1
    x(i) = (i-1)*dx;
end


% initial condition

for i = 1:m+1
  
        if x(i)>=0 && x(i)<=100
            u(1,i) = 0;
        elseif x(i)>=50 && x(i)<=220
            u(1,i) = 100*sin(pi*((x(i)-100)/120));
        else
            u(1,i) =  0;
        end
    
end

% boundary condition

for k = 1:nstep+1
    u(k,1) = 0;                                         % velocity when x = 0 at all timesteps
    u(k,m+1) = 0;                                       % velocity when x = lx at all time steps
end





t = 0;                                             % initial time
tv(1:m+1)=t;
plot3(x,tv,abs(u(1,:)),'r*')
hold on
for k = 1:nstep
    k
    t = t+dt;
  % Midpoint leapfrog method given by equation (6.18)  
   for i = 2:m
       if k == 1
        u(k+1,i)=u(k,i)+0.5*c^2*(u(k,i-1)-2*u(k,i)+u(k,i+1));   % starter solution
       else
        u(k+1,i)=2*u(k,i)-u(k-1,i)+c^2*(u(k,i-1)-2*u(k,i)+u(k,i+1));
       end
   end
    tv(1:m+1)=t;
    plot3(x,tv,abs(u(k+1,:)),'r*')
    hold on
end     
   xlabel('x')
   ylabel('t')
   zlabel('u')
% creating file and writing solution on it 

fid = fopen('CH6_P2_1_MLFM_Solution.txt','wt');
fprintf(fid,'        X     t=0.0     t=0.333   t=0.666   t=1.000  t=1.333   t=1.666   t=1.999\n\n');

for i = 1:m+1
    A = [x(i);u(1,i);u(21,i);u(41,i);u(61,i);u(81,i);u(101,i);u(121,i)];
    fprintf(fid,'%10.1f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n',A);
end
        
figure(2)
    
for i = 1:nstep+1
 
        h = plot(x(:),abs(u(i,:)),'*-r');
        set(gca, 'fontsize',15)
        xlabel('x','fontsize', 15)
        ylabel('u','fontsize', 15)
        set(h, 'XData',x(:),'YData',u(i,:));
        drawnow 
%         pause(0.1)
   
end


        
        
        
        
        
        
        
        
        
        
        
        
        
        
