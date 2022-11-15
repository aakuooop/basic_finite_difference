clc
clear all
close all

lx = 4.00;                                  % length of propagation
dx = 0.1;                                   % grid size
dt = 0.1;                                   % step size
m = lx/dx;                                  % No.of intervals
nstep = 18;                                 % Number of time step
c = dt/dx;                                  % Courant number
epsilon = 0.1;                              % damping factor
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
      
    for i = 1:m+1
        A(i) = u(k,i);
        E(i) = 0.5*u(k,i)^2;
    end
    
     a = 0;
       for i = 2:m
           a = a+1;
            d1(a) = -0.25*c*A(i-1);
            d2(a) = 1;
            d3(a) = 0.25*c*A(i+1);
        end
     % create tridiagonal matrix 
        P = diag(d1(1:m-2),-1)+diag(d2(1:m-1))+diag(d3(1:m-2),1);
        a = 0;
        for i = 2:m
            a = a+1;
    % RHS of Beam and Warming implicit method         
            if i == 2
                Q(a) = u(k,i)-0.5*c*(E(i+1)-E(i-1))...
                    +0.25*c*A(i+1)*u(k,i+1)-0.25*c*A(i-1)*u(k,i-1)...
                    +0.25*c*A(1)*u(k+1,1);
            elseif i == m
                Q(a) = u(k,i)-0.5*c*(E(i+1)-E(i-1))...
                    +0.25*c*A(i+1)*u(k,i+1)-0.25*c*A(i-1)*u(k,i-1)...
                     -0.25*c*A(m+1)*u(k+1,m+1);
            else
                Q(a) = u(k,i)-0.5*c*(E(i+1)-E(i-1))...
                    +0.25*c*A(i+1)*u(k,i+1)-0.25*c*A(i-1)*u(k,i-1)...
                    -epsilon*(u(k,i+2)-4*u(k,i+1)+6*u(k,i)...    % fourth order damping term 
                    -4*u(k,i-1)+u(k,i-2));
            end
        end
        R = P\Q';
        a = 0;
        for i = 2:m
            a = a+1;
            u(k+1,i) = R(a);
        end
        
    end
    
    fid = fopen('CH6_P9_1_The Beam and Warming Implicit Method with damping SOLUTION.txt','wt');
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
    
    legend('t=0.0','t=0.6','t=1.2','t=1.8','Location','NorthEast')
    xlabel('x')
    ylabel('u')
    
    
    
    
