clc
clear all
close all

lx = 4.00;                                  % length of propagation
dx = 0.1;                                   % grid size
dt = [0.2 0.1 0.05] ;                       % step size
m = lx/dx;                                  % No.of intervals
epsilon = 0.1;

U1 = 1;
U2 = 0;
% Grid Generation

for i = 1:m+1
    x(i) = (i-1)*dx;
end

for k = dt
    c = k/dx;
    % initial condition
    
    for i = 1:m+1
        if x(i)>=0 && x(i)<=2
            u(i) = 1;
        else
            u(i) = 0;
        end
    end
    % boundary condition
    u(1) = U1;
    u(m+1) = U2;
    
  % Descretized equation according to beam and warming method with damping
  % factor
    
    
    t = 0;
    count = 1;
    while count>0
        t = t+k
        count = 0;
        A = u ;
        E = 0.5*u.^2;
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
                Q(a) = u(i)-0.5*c*(E(i+1)-E(i-1))+0.25*c*A(i+1)*u(i+1)-0.25*c*A(i-1)*u(i-1)...
                    +0.25*c*A(1)*u(1);
            elseif i == m
                Q(a) = u(i)-0.5*c*(E(i+1)-E(i-1))+0.25*c*A(i+1)*u(i+1)-0.25*c*A(i-1)*u(i-1)...
                     -0.25*c*A(m+1)*u(m+1);
            else
                Q(a) = u(i)-0.5*c*(E(i+1)-E(i-1))+0.25*c*A(i+1)*u(i+1)-0.25*c*A(i-1)*u(i-1)...
                    -epsilon*(u(i+2)-4*u(i+1)+6*u(i)-4*u(i-1)+u(i-2)); % fourth order damping term 
            end
        end
        R = P\Q';
         
        a = 0;
        for i = 2:m
            a = a+1;
            u(i) = R(a);
        end
      
        count = count+1;
        
        if t >= 1.8
            break
        end
        
        
    end
    if k == dt(1)
        plot(x,u,'d-k')
    elseif k == dt(2)
        plot(x,u,'s-b')
    else
        plot(x,u,'*-r')
    end
    hold on
    
end

legend('c=2.0','c=1.0','c=0.5','Location','West')
xlabel('x')
ylabel('u')
