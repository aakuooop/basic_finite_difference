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
nstep = 32;                                       % No. of time steps
U1 = 0;
U2 = 0;
% Grid generation

for i = 1:m+1
    x(i) = (i-1)*dx;
end


% initial condition

for i = 1:m+1
  
        if x(i)>=0 && x(i)<=50
            u(1,i) = 0;
        elseif x(i)>=50 && x(i)<=110
            u(1,i) = 100*sin(pi*((x(i)-50)/60));
        else
            u(1,i) =  0;
        end
    
end

% boundary condition

for k = 1:nstep+1
    u(k,1) = 0;                                         % velocity when x = 0 at all timesteps
    u(k,m+1) = 0;                                       % velocity when x = lx at all time steps
end



% first upwind differencing technique according to equation (6.16)

t = 0;                                             % initial time
tv(1:m+1)=t;
plot3(x,tv,u(1,:),'r*')
hold on
for k = 1:nstep
    k
    t = t+dt
   for i = 2:m
        u(k+1,i) = u(k,i)-c*(u(k,i)-u(k,i-1));
   end
    tv(1:m+1)=t;
    plot3(x,tv,u(k+1,:),'r*')
    hold on
end     
   xlabel('x')
   ylabel('t')
   zlabel('u')
% creating file and writing solution on it 

fid = fopen('CH6_P1_1_FUDM_Solution.txt','wt');
fprintf(fid,'        X      t=0.0     t=0.10    t=0.18    t=0.28    t=0.37    t=0.45    t=0.53\n\n');

for i = 1:m+1
    A = [x(i);u(1,i);u(7,i);u(12,i);u(18,i);u(23,i);u(28,i);u(33,i)];
    fprintf(fid,'%10.1f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n',A);
end
        
figure(2)
    
for i = 1:nstep+1
 
        h = plot(x(:),u(i,:),'*-r');
        set(gca, 'fontsize',15)
        xlabel('x','fontsize', 15)
        ylabel('u','fontsize', 15)
        set(h, 'XData',x(:),'YData',u(i,:));
        drawnow 
   
end


        
        
        
        
        
        
        
        
        
        
        
        
        
        
