clc
clear all
close all

% Given input data
a = 200;                                          % speed of sound
dx =1.0;                                          % Grid size
dt = 0.0025;                                      % step size
c = a*dt/dx;                                      % Courant number
lx = 70.00;                                       % length of one-dimensional tube
m = lx/dx;                                        % No. of intervals
nstep = 60;                                       % No. of time steps
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


% Second-order accurate Lax-Wendroff method given by equation (6.5)

t = 0;                                                    % initial time
tv(1:m+1)=t;
plot3(x,tv,u(1,:),'r*')
hold on
for k = 1:nstep
    k
    t = t+dt
   for i = 2:m
        u(k+1,i) = u(k,i)-c*0.5*(u(k,i+1)-u(k,i-1))+0.5*c^2*(u(k,i+1)-2*u(k,i)+u(k,i-1));
   end
    tv(1:m+1)=t;
plot3(x,tv,u(k+1,:),'r*')
hold on
end
xlabel('x')
ylabel('t')
zlabel('u')
   
% % creating file and writing solution on it 

fid = fopen('P6_1_2_LWM_Solution.txt','wt');
fprintf(fid,'        X      t=0.0     t=0.025   t=0.05   t=0.075    t=0.1    t=0.125    t=0.15\n\n');

for i = 1:m+1
    A = [x(i);u(1,i);u(11,i);u(21,i);u(31,i);u(41,i);u(51,i);u(61,i)];
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


        
        
        
        
        
        
        
        
        
        
        
        
        
        
