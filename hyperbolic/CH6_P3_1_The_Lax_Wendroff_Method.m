clc
clear all
close all

lx = 4.00;                                  % length of propagation
dx = 0.1;                                   % grid size
dt = 0.1;                                   % step size
m = lx/dx;                                  % No.of intervals 
nstep = 18;                                 % Number of time step 
c = dt/dx;

% Grid Generation 

for i = 1:m+1
    x(i) = (i-1)*dx;
end


% initial condition 

for i = 1:m+1
    if x(i)>=0 && x(i)<=2
        u(1,i) = 1;
    else
        u(1,i) = 0;
    end
end


% boundary condition 

for k = 1:nstep+1
    u(k,1) = 1;
    u(k,m+1) = 0;
end


t = 0;
for k = 1:nstep
    k
    t = t+dt
    % Descretized equation according to Lax-Wendroff  method 
        for i = 2:m
            u(k+1,i)=u(k,i)-0.25*c*(u(k,i+1)^2-u(k,i-1)^2)+0.125*c^2*((u(k,i+1)+u(k,i))*(u(k,i+1)^2-u(k,i)^2)-(u(k,i)...
                +u(k,i-1))*(u(k,i)^2-u(k,i-1)^2));
        end
end

fid = fopen('CH6_P3_1_THE LAX WENDROFF METHOD_SOLUTION.txt','wt');
fprintf(fid,'       x      t=0.0     t=0.3     t=0.6     t=0.9     t=1.2     t=1.5     t=1.8\n\n');

for i = 1:m+1
    A = [x(i);u(1,i);u(4,i);u(7,i);u(10,i);u(13,i);u(16,i);u(19,i)];
    fprintf(fid,'%10.2f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',A);
end
plot(x,u(1,:),'r*-')
hold on
plot(x,u(7,:),'bs-')
hold on
plot(x,u(13,:),'kd-')
hold on
plot(x,u(19,:),'g+-')
hold off

legend('t=0.0','t=0.6','t=1.2','t=1.8','Location','West')
xlabel('x')
ylabel('u')




