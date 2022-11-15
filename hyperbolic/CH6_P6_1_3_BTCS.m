clc
clear all
close all

% Given input data


a = 200;                                          % speed of sound 
dx =1.0;                                          % Grid size
dt = 0.00125;                                     % step size
c = a*dt/dx;                                      % Courant number
lx = 70.00;                                       % length of one-dimensional tube
m = lx/dx;                                        % No. of intervals
nstep = 120;                                      % No. of time steps
U1 = 0;
U2 = 0;
% Grid generation

for i = 1:m+1
    x(i) = (i-1)*dx;
end


% initial condition
theta = atan(2);
for i = 1:m+1
  
        if x(i)>=0 && x(i)<=5
            u(1,i) = 0;
        elseif x(i)>5 && x(i)<=15
            u(1,i) = (x(i)-5)*tan(theta);
        elseif x(i)>15 && x(i)<=25
            u(1,i) = (25-x(i))*tan(theta);
        else
            u(1,i) =  0;
        end
    
end
% boundary condition

for k = 1:nstep+1
    u(k,1) = 0;                                         % velocity when x = 0 at all timesteps
    u(k,m+1) = 0;                                       % velocity when x = lx at all time steps
end


%  Euler's BTCS method given by equation (6.6)

t = 0;                                             % initial time
tv(1:m+1)=t;
plot3(x,tv,u(1,:),'r*')
hold on
xlabel('x')
ylabel('t')
zlabel('u')

for i = 1:m-1
    d1(i)=0.5*c;
    d2(i)=-1;
    d3(i)=-0.5*c;
end

% create tridiagonal matrix 
P = diag(d1(1:m-2),-1)+diag(d2(1:m-1))+diag(d3(1:m-2),1);
for k = 1:nstep
    t = t+dt;
    r = 0;
   for i = 2:m
       r = r+1;
       
   % RHS of equation 6.6     
       if i ==2
       Q(r) = -u(k,i)-0.5*c*u(k+1,i-1);
       elseif i==m
           Q(r) = -u(k,i)+0.5*c*u(k+1,i+1);
       else
            Q(r) = -u(k,i);
       end
   end
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

fid = fopen('P6_1_3_BTCS_Solution.txt','wt');
fprintf(fid,'        X      t=0.0     t=0.025   t=0.05   t=0.075    t=0.1    t=0.125    t=0.15\n\n');

for i = 1:m+1
    A = [x(i);u(1,i);u(21,i);u(41,i);u(61,i);u(81,i);u(101,i);u(121,i)];
    fprintf(fid,'%10.1f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n',A);
end
        
figure(2)
    
for k = 1:nstep+1
 
        h = plot(x,u(k,:),'*-r');
        set(gca, 'fontsize',15)
        xlabel('x','fontsize', 15)
        ylabel('u','fontsize', 15)
        set(h, 'XData',x(:),'YData',u(k,:));
        drawnow 
   
end


        
        
        
        
        
        
        
        
        
        
        
        
        
        
