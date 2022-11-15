clc
clear all
close all

% Given input data


a = 300;                                          % speed of sound at the altitude of 10 km
dx =5.0;                                          % Grid size
dt = 0.015;                                       % step size
c = a*dt/dx;                                      % Courant number
lx = 300;                                         % length of one-dimensional tube
m = lx/dx;                                        % No. of intervals
nstep = 36;                                       % No. of time steps
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



t = 0;                                             % initial time
tv(1:m+1)=t;
plot3(x,tv,u(1,:),'r*')
hold on

% Fully implicit Euler's BTCS method given by equation (6.6)

for i = 1:m-1
    d1(i)=0.5*c;
    d2(i)=-1;
    d3(i)=-0.5*c;
end

% Formulate tridiagonal matrix 
P = diag(d1(1:m-2),-1)+diag(d2(1:m-1))+diag(d3(1:m-2),1);
for k = 1:nstep
    t = t+dt;
    r = 0;
   for i = 2:m
       r = r+1;
       if i ==2
       Q(r) = -u(k,i)-0.5*c*u(k+1,i-1);
       elseif i==m
           Q(r) = -u(k,i)+0.5*c*u(k+1,i+1);
       else
            Q(r) = -u(k,i);
       end
   end
 % solving matrix
    R = P\Q';
    a = 0;
    for i = 2:m
        a = a+1;
        u(k+1,i)= R(a);
    end
    tv(1:m+1)=t;
plot3(x,tv,u(k+1,:),'r*')
hold on
    
end    
xlabel('x')
ylabel('t')
zlabel('u')
   
% creating file and writing solution on it 

fid = fopen('CH6_P1_1_BTCS_Solution.txt','wt');
fprintf(fid,'        X      t=0.0     t=0.09    t=0.18    t=0.27    t=0.36    t=0.45    t=0.54\n\n');

for i = 1:m+1
    A = [x(i);u(1,i);u(7,i);u(13,i);u(19,i);u(25,i);u(31,i);u(37,i)];
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
        
        
        
        
        
        
        
        
        
        
