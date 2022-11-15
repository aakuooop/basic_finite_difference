clc
clear all
close all

lx = 40.00;                                  % length of propagation
dx = 1;                                   % grid size
dt = 0.2;                                   % step size
m = lx/dx;                                  % No.of intervals 
nstep = 12;                                 % Number of time step 
c = dt/dx;

% Grid Generation 

for i = 1:m+1
    x(i) = (i-1)*dx;
end


% initial condition 

for i = 1:m+1
    if x(i)>=0 && x(i)<=20
        u(1,i) = 5.0;
    else
        u(1,i) = 0;
    end
end


% boundary condition 

for k = 1:nstep+1
    u(k,1) = 5;
    u(k,m+1) = 0;
end

% Descretized equation according to Lax-Wendroff  method 
t = 0;
for k = 1:nstep
    k
    t = t+dt
        for i = 2:m
            u(k+1,i)=u(k,i)-0.25*c*(u(k,i+1)^2-u(k,i-1)^2)+0.125*c^2*((u(k,i+1)+u(k,i))*(u(k,i+1)^2-u(k,i)^2)-(u(k,i)...
                +u(k,i-1))*(u(k,i)^2-u(k,i-1)^2));
        end
end

fid = fopen('P6_2_1_LAX-WENDROFF METHOD_SOLUTION.txt','wt');
fprintf(fid,'       x      t=0.0     t=0.4     t=0.8     t=1.2     t=1.6     t=2.0     t=2.4\n\n');

for i = 1:m+1
    A = [x(i);u(1,i);u(3,i);u(5,i);u(7,i);u(9,i);u(11,i);u(13,i)];
    fprintf(fid,'%10.2f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',A);
end
plot(x,u(1,:),'r*-')
hold on
plot(x,u(5,:),'b^-')
hold on
plot(x,u(9,:),'kd-')
hold on
plot(x,u(13,:),'go-')
hold off

legend('t=0.0','t=0.8','t=1.6','t=2.4','Location','West')
xlabel('x')
ylabel('u')




