clc
clear all
close all

lx = 4.00;                                  % length of propagation
dx = 0.05;                                   % grid size
dt = 0.01;                                   % step size
m = lx/dx;                                  % No.of intervals 
nstep = 600 ;                                 % Number of time step 
c = dt/dx;

% Grid Generation 

for i = 1:m+1
    x(i) = (i-1)*dx;
end


% initial condition 

for i = 1:m+1
    if x(i)>=0 && x(i)<0.25
        u(1,i) = 1.0;
    elseif x(i)>=0.25 && x(i)<=1.25
        u(1,i) = 1.25-x(i)
    else
        u(1,i) = 0;
    end
end


% boundary condition 

for k = 1:nstep+1
    u(k,1) = 1.0;
    u(k,m+1) = 0.0;
end

% Descretized equation according to Lax-Wendroff  method 
t = 0;
for k = 1:nstep
    k
    t = t+dt
        for i = 2:m
            u(k+1,i)=0.5*(u(k,i+1)+u(k,i-1))-0.25*c*((u(k,i+1))^2-(u(k,i-1))^2);
        end
end

fid = fopen('P6_3_1_THE Lax-Wendroff SOLUTION.txt','wt');
fprintf(fid,'       x      t=0.0     t=2.0     t=4.0     t=6.0\n\n');

for i = 1:2:m+1
    A = [x(i);u(1,i);u(201,i);u(401,i);u(601,i)];
    fprintf(fid,'%10.2f%10.5f%10.5f%10.5f%10.5f\n',A);
end
plot(x,u(1,:),'r*-')
hold on
plot(x,u(201,:),'b^-')
hold on
plot(x,u(401,:),'kd-')
hold on
plot(x,u(601,:),'go-')
hold off

legend('t=0.0','t=2.0','t=4.0','t=6.0','Location','West')
xlabel('x')
ylabel('u')

