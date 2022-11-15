clc
clear all
close all

lx = 40.00;                                  % length of propagation
dx = 1;                                      % grid size
dt = 0.1;                                   % step size
m = lx/dx;                                  % No.of intervals 
nstep = 24;                                 % Number of time step 
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
    u(k,1) = 5.0;
    u(k,m+1) = 0;
end
u_predict(1) = 5.0;
u_predict(m+1) = 0;

% Descretized equation according to MacCormack method 
t = 0;
for k = 1:nstep
    k
    t = t+dt
    
    % predictor 
       for i = 2:m
           u_predict(i) =  u(k,i)-0.5*c*(u(k,i+1)^2-u(k,i)^2);
       end
    % corrector
        for i = 2:m
           u(k+1,i) = 0.5*(u(k,i)+u_predict(i)-0.5*c*(u_predict(i)^2-u_predict(i-1)^2));
        end
end

fid = fopen('P6_2_1_MacCormack.txt','wt');
fprintf(fid,'       x      t=0.0     t=0.4     t=0.8     t=1.2     t=1.6     t=2.0     t=2.4\n\n');

for i = 1:2:m+1
    A = [x(i);u(1,i);u(5,i);u(9,i);u(13,i);u(17,i);u(21,i);u(25,i)];
    fprintf(fid,'%10.2f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',A);
end
plot(x,u(1,:),'r*-')
hold on
plot(x,u(9,:),'b^-')
hold on
plot(x,u(17,:),'kd-')
hold on
plot(x,u(25,:),'go-')
hold off

legend('t=0.0','t=0.8','t=1.6','t=2.4','Location','West')
xlabel('x')
ylabel('u')
